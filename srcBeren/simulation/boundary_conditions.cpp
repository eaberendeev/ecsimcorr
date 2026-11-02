#include "simulation/boundary_conditions.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "containers.h"
#include "grid/domain.h"
#include "particle_array.h"
#include "timer.h"

namespace {
// Кандидат на вторичную эмиссию: частица, поглощённая на грани face_, и ключ
// канонического порядка (линейный индекс ячейки), чтобы эмиссия выполнялась
// детерминированно независимо от числа потоков и расписания.
struct EmissionCandidate {
    Particle p;
    Face face;
    int64_t order_key;
};
}   // namespace

ParticleFateResult BoundaryConditionHandler::classify_particle(const Particle& p, const ParticlesArray& particles,
                                                               const Domain& domain, Face hit_face) const {
    for (const ParticleBoundaryCondition* condition : particle_conditions_) {
        const ParticleFateResult result = condition->classify_particle(p, particles, domain, hit_face);
        if (result.fate != ParticleFate::Unhandled)
            return result;
    }
    return {};
}

void BoundaryConditionHandler::apply_to_particles(
    ParticlesArray& particles, std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species,
    const Domain& domain, ParticleBoundaryEvents* captured) {
    RECORD_TIMER;

    Array3D<std::vector<Particle>>& data = particles.cells();
    int nx = data.size().x(), ny = data.size().y(), nz = data.size().z();

    // Двухфазный подход без #pragma omp critical:
    // 1) параллельная фаза — условия классифицируют частицы и пишут только в
    //    локальные (поточные) буферы эмиттера, диагностики и кандидатов;
    // 2) последовательная фаза — слияние поточных буферов и детерминированная
    //    (по каноническому порядку ячеек) вторичная эмиссия.
    const int nthreads = omp_get_max_threads();
    std::vector<BoundaryEmitter> local_emitters(nthreads);
    std::vector<SpeciesStepState> local_diag(nthreads);
    std::vector<std::vector<EmissionCandidate>> local_candidates(nthreads);

    // Сортируем один раз до параллельного поиска: записи завершены после
    // барьера депозита, поэтому бинарный поиск не гоняется с ними.
    if (captured)
        captured->sort_for_lookup();

#pragma omp parallel
    {
        const int tid = omp_get_thread_num();

        // Классификатор частицы: возвращает true, если частица должна быть
        // удалена из ячейки. Не содержит никакого разделяемого состояния — все
        // записи идут в поточные буферы local_*[tid].
        const auto classifyAndRecord = [&](const Particle& p, int ix, int iy, int iz) -> bool {
            if (domain.contains(p.coord))
                return false;

            if (periodic_[0] || periodic_[1] || periodic_[2]) {
                Particle p_new = p;
                p_new.coord = wrap_periodic(p.coord, domain);

                if (domain.contains(p_new.coord)) {
                    local_emitters[tid].emit_current_species(p_new);
                    return true;
                }
            }

            const ParticleFateResult* saved = captured ? captured->find(p) : nullptr;
            const ParticleFateResult res = saved ? *saved : classify_particle(p, particles, domain);
            if (res.fate == ParticleFate::Removed) {
                local_diag[tid].add_loss(res.face, get_energy_particle(p.velocity, particles.mass(), particles.mpw()));
                if (!emissions_[static_cast<size_t>(res.face)].empty())
                    local_candidates[tid].push_back({p, res.face, (static_cast<int64_t>(ix) * ny + iy) * nz + iz});
            } else if (res.fate == ParticleFate::Reflected) {
                Particle reflected = p;
                reflected.coord = domain.geom.reflect_from_face(res.face, p.coord);
                reflected.velocity.z() = -p.velocity.z();
                local_emitters[tid].emit_current_species(reflected);
                local_diag[tid].add_reflected(res.face);
            }
            return true;
        };

#pragma omp for collapse(3) schedule(dynamic, 32)
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    std::vector<Particle>& cell = data(ix, iy, iz);

                    // Walk backwards: a swapped survivor was already processed.
                    // Its saved event need not move with the particle's storage address.
                    for (int i = static_cast<int>(cell.size()) - 1; i >= 0; --i) {
                        if (classifyAndRecord(cell[i], ix, iy, iz)) {
                            std::swap(cell[i], cell.back());
                            cell.pop_back();
                        }
                    }
                }
            }
        }
    }

    // Последовательная фаза.
    timer::commonTimer timerEmitter("add particles");
    // добавить новые частицы текущего сорта (отражённые и замкнутые по периодике)
    for (int t = 0; t < nthreads; ++t) {
        for (const Particle& p : local_emitters[t].current_species_particles()) {
            particles.add_particle(p);
        }
        // слияние поточной диагностики в общую
        particles.step_state.merge_from(local_diag[t]);
        // Страховка: в параллельной фазе условиям разрешена только эмиссия в
        // ТЕКУЩИЙ сорт (вызовы emit_to_species потерялись бы). Если условие
        // всё же записало в species_buffers_ (нарушение контракта), не
        // теряем частицы — переносим в общий эмиттер; flush_species в конце
        // шага разложит их по сортам.
        for (const auto& [name, buf] : local_emitters[t].other_species_particles()) {
            for (const Particle& p : buf) {
                emitter.emit_to_species(name, p);
            }
        }
    }

    // Вторичная эмиссия: собираем всех кандидатов и сортируем по каноническому
    // порядку ячеек, поэтому порядок эмиссии детерминирован независимо от числа
    // потоков и расписания.
    std::vector<EmissionCandidate> candidates;
    for (int t = 0; t < nthreads; ++t) {
        candidates.insert(candidates.end(), local_candidates[t].begin(), local_candidates[t].end());
    }
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const EmissionCandidate& a, const EmissionCandidate& b) { return a.order_key < b.order_key; });
    for (const EmissionCandidate& c : candidates) {
        for (SecondaryEmissionModel& m : emissions_[static_cast<size_t>(c.face)]) {
            m.emit(c.p, particles, all_species, emitter, domain);
        }
    }
}

void BoundaryConditionHandler::flush_species(
    std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) {
    RECORD_TIMER;
    // добавить частицы в другие сорта по имени
    const auto& other = emitter.other_species_particles();

    for (const auto& kv : other) {
        const std::string& name = kv.first;
        if (all_species.count(name) == 0) {
            throw std::runtime_error("Species " + name + " not found when adding from emitter");
        }
        auto& target_species = all_species[name];
        for (const Particle& p : kv.second) {
            target_species->add_particle(p);
        }
    }
    emitter.clear_other_species_buffers();
}

void BoundaryConditionHandler::validate_emissions(
    const std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) const {
    RECORD_TIMER;
    // Проверяем, что все сорта-продукты вторичной эмиссии существуют среди
    // зарегистрированных сортов частиц. Вызывается ПОСЛЕ инициализации всех
    // сортов и до основного цикла — бросать здесь безопасно (последовательный
    // код).
    std::string missing;
    for (size_t f = 0; f < emissions_.size(); ++f) {
        const auto face = static_cast<Face>(f);
        for (const SecondaryEmissionModel& m : emissions_[f]) {
            if (!all_species.count(m.product())) {
                missing +=
                    "  face " + std::to_string(static_cast<int>(face)) + ": product species \"" + m.product() + "\"\n";
            }
        }
    }
    if (!missing.empty()) {
        throw std::runtime_error(
            "\"second_emission\" references product species not defined in the particle configuration:\n" + missing);
    }
}
