#include <omp.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "boundary_conditions.h"
#include "containers.h"
#include "timer.h"

namespace {
// Кандидат на вторичную эмиссию: частица, поглощённая на грани face_, и ключ
// канонического порядка (линейный индекс ячейки), чтобы эмиссия выполнялась
// детерминированно независимо от числа потоков и расписания.
struct EmissionCandidate {
    Particle p;
    Face face;
    long long order_key;
};
}   // namespace

void BoundaryConditionHandler::apply_to_particles(
    ParticlesArray& particles, std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species,
    const Domain& domain) {
    RECORD_TIMER;

    auto& data = particles.particlesData;
    int nx = data.size().x(), ny = data.size().y(), nz = data.size().z();

    // Двухфазный подход без #pragma omp critical:
    // 1) параллельная фаза — условия классифицируют частицы и пишут только в
    //    локальные (поточные) буферы эмиттера, диагностики и кандидатов;
    // 2) последовательная фаза — слияние поточных буферов и детерминированная
    //    (по каноническому порядку ячеек) вторичная эмиссия.
    const int nthreads = omp_get_max_threads();
    std::vector<BoundaryEmitter> local_emitters(nthreads);
    std::vector<SpeciesDiagStats> local_diag(nthreads);
    std::vector<std::vector<EmissionCandidate>> local_candidates(nthreads);

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

            for (const auto& cond : conditions_) {
                auto res = cond->apply_to_particle(p, particles, local_emitters[tid], domain);
                if (res.fate != ParticleFate::Unhandled) {
                    if (res.fate == ParticleFate::Removed) {
                        local_diag[tid].add_loss(res.face,
                                                 get_energy_particle(p.velocity, particles.mass(), particles.mpw()));
                        if (emissions_.contains(res.face)) {
                            local_candidates[tid].push_back({p, res.face, ((long long) ix * ny + iy) * nz + iz});
                        }
                    } else if (res.fate == ParticleFate::Reflected) {
                        local_diag[tid].add_reflected(res.face);
                    }
                    break;
                }
            }
            return true;
        };

#pragma omp for collapse(3) schedule(dynamic, 32)
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    std::vector<Particle>& cell = data(ix, iy, iz);

                    for (int i = 0; i < std::ssize(cell);) {
                        if (classifyAndRecord(cell[i], ix, iy, iz)) {
                            std::swap(cell[i], cell[cell.size() - 1]);
                            cell.resize(cell.size() - 1);
                        } else {
                            i += 1;
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
        for (const auto& p : local_emitters[t].current_species_particles()) {
            particles.add_particle(p);
        }
        // слияние поточной диагностики в общую
        particles.diag.merge_from(local_diag[t]);
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
    for (const auto& c : candidates) {
        auto it = emissions_.find(c.face);
        if (it == emissions_.end())
            continue;
        for (auto& m : it->second) {
            m.emit(c.p, particles, all_species, emitter, domain);
        }
    }
}

void BoundaryConditionHandler::flush_species(
    std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) {
    // добавить частицы в другие сорта по имени
    const auto& other = emitter.other_species_particles();

    for (const auto& kv : other) {
        const std::string& name = kv.first;
        if (all_species.count(name) == 0) {
            throw std::runtime_error("Species " + name + " not found when adding from emitter");
        }
        auto& target_species = all_species[name];
        for (const auto& p : kv.second) {
            target_species->add_particle(p);
        }
    }
    emitter.clear_other_species_buffers();
}

void BoundaryConditionHandler::validate_emissions(
    const std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>& all_species) const {
    // Проверяем, что все сорта-продукты вторичной эмиссии существуют среди
    // зарегистрированных сортов частиц. Вызывается ПОСЛЕ инициализации всех
    // сортов и до основного цикла — бросать здесь безопасно (последовательный
    // код).
    std::string missing;
    for (const auto& [face, models] : emissions_) {
        for (const auto& m : models) {
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

void BoundaryConditionHandler::load_from_json(const nlohmann::json& sys_config, const Domain& domain) {
    conditions_.clear();
    emissions_.clear();
    if (!sys_config.contains("Boundary_conditions"))
        return;
    const auto& config = sys_config["Boundary_conditions"];
    if (!config.is_array()) {
        std::cerr << "BoundaryConditionHandler: expected array in JSON\n";
        return;
    }

    for (const auto& item : config) {
        if (item.empty())
            continue;

        auto it = item.begin();
        std::string type = it.key();
        const auto& params = it.value();

        // "open" accepts only an array of faces; all faces are gathered into
        // a single OpenBoundaryConditionArray (this class exists exactly for
        // such batching)
        if (type == "open") {
            if (!params.is_array()) {
                throw std::runtime_error(
                    "\"open\" expects an array of faces, e.g. \"open\": [{\"face\": \"XMIN\"}, {\"face\": \"XMAX\"}]");
            }
            std::vector<Face> faces;
            for (const nlohmann::json& faceJs : params) {
                faces.push_back(string_to_face(faceJs.at("face")));
            }
            conditions_.push_back(std::make_unique<OpenBoundaryConditionArray>(faces));
            continue;
        }

        // Other types accept either a single face object or an array of them
        // (e.g. {"periodic": [{"face": "YMIN"}, {"face": "YMAX"}]})
        if (params.is_array()) {
            for (const nlohmann::json& faceJs : params) {
                add_condition(type, faceJs, domain);
            }
        } else {
            add_condition(type, params, domain);
        }
    }

    // Проверка ПОСЛЕ обработки всех записей: для каждой грани, на которой
    // задана вторичная эмиссия, должно существовать consuming-условие (open,
    // bphi и т.п.) на той же грани — иначе эмиссия никогда не сработает
    // (частицы не будут удалены на этой грани). Условие может идти в JSON
    // после second_emission, поэтому проверяем только в конце цикла.
    for (const auto& [face, models] : emissions_) {
        bool has_consumer = false;
        for (const auto& cond : conditions_) {
            if (cond->face() == face) {
                has_consumer = true;
                break;
            }
            // OpenBoundaryConditionArray покрывает сразу несколько граней.
            if (const auto* arr = dynamic_cast<const OpenBoundaryConditionArray*>(cond.get())) {
                for (int i = 0; i < arr->createdBc_; ++i) {
                    if (arr->faces_[i] == face) {
                        has_consumer = true;
                        break;
                    }
                }
            }
            if (has_consumer)
                break;
        }
        if (!has_consumer) {
            throw std::runtime_error("\"second_emission\" on face " + std::to_string(static_cast<int>(face)) +
                                     " needs a consuming boundary condition (e.g. \"open\" or \"bphi\") on the same "
                                     "face, otherwise it never fires.");
        }
    }
}

void BoundaryConditionHandler::add_condition(const std::string& type, const nlohmann::json& params,
                                             const Domain& domain) {
    const std::string face_str = params.at("face");
    const Face face = string_to_face(face_str);

    if (face == Face::CYLINDER && !domain.geom.use_cylinder) {
        throw std::runtime_error(
            "Boundary condition uses face \"CYLINDER\", but \"CylinderDomain\" is not configured.\n"
            "Add to system_config:\n"
            "  \"CylinderDomain\": {\"radius\": <value>, \"center\": [<x>, <y>]}");
    }

    if (type == "bphi") {
        const double electron_threshold_energy = params.at("electron_threshold_energy_");
        const double radius = params.at("radius");
        const double gap = params.at("gap");
        conditions_.push_back(std::make_unique<BphiCondition>(face, gap, radius, electron_threshold_energy));
    } else if (type == "second_emission") {
        const std::string product = params.value("product", "Electrons");
        if (!params.contains("sources") || !params["sources"].is_array() || params["sources"].empty()) {
            throw std::runtime_error(
                "\"second_emission\" requires a non-empty \"sources\" array, e.g.:\n"
                "  {\"second_emission\": {\"face\": \"ZMIN\", \"product\": \"Electrons\", \"sources\": [\n"
                "      {\"species\": \"Ions\", \"yield\": 0.3, \"energy\": {\"type\": \"fixed\", \"kev\": 2.0}},\n"
                "      {\"species\": \"Electrons\", \"yield\": 0.1, \"energy\": {\"type\": \"temperature\", "
                "\"temperature_kev\": [1.0,1.0,1.0]}}\n"
                "      {\"species\": \"Ions2\", \"yield\": 0.5, \"energy\": {\"type\": \"fraction\", \"fraction\": "
                "0.5}}\n"
                "    ]}}");
        }
        std::vector<EmissionSourceRule> rules;
        for (const auto& src : params["sources"]) {
            EmissionSourceRule rule;
            rule.species = src.at("species").get<std::string>();
            if (src.contains("yield") && src["yield"].is_number()) {
                rule.yield_model = EmissionSourceRule::YieldModel::Constant;
                rule.yield = src["yield"].get<double>();
                if (rule.yield < 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"yield\" must be >= 0, got " + std::to_string(rule.yield));
                }
            } else if (src.contains("yield") && src["yield"].is_object()) {
                const nlohmann::json& ycfg = src["yield"];
                const std::string ytype = ycfg.value("type", "");
                if (ytype == "vaughan") {
                    rule.yield_model = EmissionSourceRule::YieldModel::Vaughan;
                    rule.delta_max = ycfg.at("delta_max").get<double>();
                    rule.energy_max_kev = ycfg.at("energy_max_kev").get<double>();
                    rule.threshold_kev = ycfg.value("threshold_kev", 0.0);
                    if (rule.delta_max <= 0.0) {
                        throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                                 "\": \"delta_max\" must be > 0, got " +
                                                 std::to_string(rule.delta_max));
                    }
                    if (rule.energy_max_kev <= 0.0) {
                        throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                                 "\": \"energy_max_kev\" must be > 0, got " +
                                                 std::to_string(rule.energy_max_kev));
                    }
                    if (rule.threshold_kev < 0.0) {
                        throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                                 "\": \"threshold_kev\" must be >= 0, got " +
                                                 std::to_string(rule.threshold_kev));
                    }
                } else if (ytype == "threshold") {
                    rule.yield_model = EmissionSourceRule::YieldModel::Threshold;
                    rule.yield = ycfg.at("yield").get<double>();
                    rule.threshold_kev = ycfg.at("threshold_kev").get<double>();
                    if (rule.yield < 0.0) {
                        throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                                 "\": \"yield\" must be >= 0, got " + std::to_string(rule.yield));
                    }
                    if (rule.threshold_kev < 0.0) {
                        throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                                 "\": \"threshold_kev\" must be >= 0, got " +
                                                 std::to_string(rule.threshold_kev));
                    }
                } else {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": unknown \"yield\" type \"" + ytype +
                                             "\" (expected number, {\"type\":\"vaughan\",...} or "
                                             "{\"type\":\"threshold\",...})");
                }
            } else {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                         "\": \"yield\" must be a number or an object");
            }
            if (!src.contains("energy") || !src["energy"].is_object()) {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                         "\" requires an \"energy\" object, e.g. {\"type\": \"fixed\", \"kev\": 2.0}");
            }
            const auto& energy_cfg = src["energy"];
            const std::string etype = energy_cfg.at("type").get<std::string>();
            if (etype == "fixed") {
                rule.energy_type = EmissionSourceRule::EnergyType::Fixed;
                rule.fixed_kev = energy_cfg.at("kev").get<double>();
                if (rule.fixed_kev < 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"kev\" must be >= 0, got " + std::to_string(rule.fixed_kev));
                }
            } else if (etype == "temperature") {
                rule.energy_type = EmissionSourceRule::EnergyType::Temperature;
                rule.temperature_kev = util::parse_double3(energy_cfg.at("temperature_kev"));
                if (rule.temperature_kev.x() < 0.0 || rule.temperature_kev.y() < 0.0 ||
                    rule.temperature_kev.z() < 0.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": every component of \"temperature_kev\" must be >= 0, got [" +
                                             std::to_string(rule.temperature_kev.x()) + ", " +
                                             std::to_string(rule.temperature_kev.y()) + ", " +
                                             std::to_string(rule.temperature_kev.z()) + "]");
                }
                rule.temperature_mean =
                    energy_cfg.contains("mean") ? util::parse_double3(energy_cfg["mean"]) : Vector3R(0, 0, 0);
            } else if (etype == "fraction") {
                rule.energy_type = EmissionSourceRule::EnergyType::Fraction;
                rule.fraction = energy_cfg.at("fraction").get<double>();
                if (rule.fraction < 0.0 || rule.fraction > 1.0) {
                    throw std::runtime_error("\"second_emission\" source \"" + rule.species +
                                             "\": \"fraction\" must be in [0..1], got " +
                                             std::to_string(rule.fraction));
                }
            } else {
                throw std::runtime_error("\"second_emission\" source \"" + rule.species + "\": unknown energy type \"" +
                                         etype + "\" (expected \"fixed\", \"temperature\" or \"fraction\")");
            }
            rules.push_back(rule);
        }
        for (const auto& rule : rules) {
            if (rule.yield_model == EmissionSourceRule::YieldModel::Constant && rule.species == product &&
                rule.yield >= 1.0) {
                throw std::runtime_error("\"second_emission\": source \"" + rule.species + "\" equals product \"" +
                                         product + "\" with constant yield >= 1 — self-amplifying emission would run away; use \"vaughan\"/\"threshold\" with an energy gate or yield < 1");
            }
        }
        // ВАЖНО: эмиссия — НЕ consuming-условие: ничего не добавляем в
        // conditions_, только в emissions_.
        emissions_[face].emplace_back(face, product, std::move(rules));
    } else if (type == "second_emisson") {
        // Устаревший вариант написания (опечатка): сообщаем и просим переименовать.
        throw std::runtime_error(
            "Boundary condition \"second_emisson\" is a typo and no longer supported.\n"
            "Rename it to \"second_emission\" and use the new format, e.g.:\n"
            "  {\"second_emission\": {\"face\": \"ZMIN\", \"product\": \"Electrons\", \"sources\": [\n"
            "      {\"species\": \"Ions\", \"yield\": 0.3, \"energy\": {\"type\": \"fixed\", \"kev\": 2.0}}]}}");
    } else if (type == "periodic") {
        // Periodic folding of fields/operators acts on the whole axis at once
        // (it pairs the low and high boundary layers mutually), so a single
        // condition per axis is sufficient. XMIN+XMAX entries in the config
        // describe the same axis; adding both would apply the fold twice and
        // double the field values in the boundary layers.
        int axis = -1;
        if (face == Face::XMIN || face == Face::XMAX) {
            axis = 0;
        } else if (face == Face::YMIN || face == Face::YMAX) {
            axis = 1;
        } else if (face == Face::ZMIN || face == Face::ZMAX) {
            axis = 2;
        }
        if (axis >= 0 && periodic_[axis]) {
            std::cout << "BoundaryConditionHandler: duplicate periodic condition on face " << face_str
                      << " is ignored (axis already periodic)\n";
            return;
        }
        conditions_.push_back(std::make_unique<PeriodicBoundaryCondition>(face));
        if (axis >= 0) {
            periodic_[axis] = true;
        }
    } else if (type == "electron_reflection") {
        const double radius = params.value("radius", 5.0);
        const double energy_threshold = params.value("energy_threshold", 0.1) / SGS::MC2;
        conditions_.push_back(std::make_unique<ElectronReflectionCondition>(face, radius, energy_threshold));
    } else if (type == "er0") {
        const double inner_radius = params.value("inner_radius", 5.0);
        const double width = params.value("width", 2.0);
        const double potential_drop = params.value("potential_drop", 0.1) / SGS::MC2;
        conditions_.push_back(std::make_unique<Er0Condition>(face, inner_radius, width, potential_drop));
    } else {
        std::cerr << "Unknown boundary condition type: " << type << std::endl;
    }
}
