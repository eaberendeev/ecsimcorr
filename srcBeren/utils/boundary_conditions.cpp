#include "boundary_conditions.h"

#include <functional>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "ParticlesArray.h"
#include "World.h"
#include "containers.h"

// -----------------------------------------------
// Конкретные реализации
// -----------------------------------------------

void PeriodicBoundaryCondition::apply_to_particle(
    const Particle& p, const std::string& /*species_name*/,
    BoundaryEmitter& emitter, const Domain& domain) {
    Particle p_new = p;
    //domain.make_point_periodic(p_new.coord);
    emitter.emit_current_species(p_new);
}

void SecondEmissionCondition::apply_to_particle(const Particle& p,
                                                const std::string& species_name,
                                                BoundaryEmitter& emitter,
                                                const Domain& domain) {
    // Условие: вторичная эмиссия работает только для электронов
    if (species_name != "Ions") {
        return;
    }

    if (!domain.geom.is_outside_face(face_, p.coord) && is_outside_other_faces(p.coord, domain)) return;

    Particle p_new = p;
    p_new.velocity = gauss_.sample(pulse_gen_);

    if(face_ == Face::ZMIN){
        p_new.coord.z() = -p_new.coord.z();

        p_new.velocity.z() = std::abs(p_new.velocity.z());
    } else if (face_ == Face::ZMAX) {
        p_new.coord.z() = 2 * domain.cell_size().z() * domain.num_cells().z() -
                          p_new.coord.z();
        p_new.velocity.z() = std::abs(p_new.velocity.z());
    } else {
        return;
    }

    emitter.emit_to_species("Electrons", p_new);
}

void BoundaryConditionHandler::apply_to_particles(
    ParticlesArray& particles,
    [[maybe_unused]] std::unordered_map<
        std::string, std::unique_ptr<ParticlesArray>>& all_species,
    const Domain& domain) {
    double time0 = omp_get_wtime();
    auto& data = particles.particlesData;
    int nx = data.size().x(), ny = data.size().y(), nz = data.size().z();


    // имя текущего сорта
    const std::string& species_name = particles.name();
    std::cout << "Applying boundary conditions to " << species_name << "\n";
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                auto& cell = data(ix, iy, iz);

                auto it = std::remove_if(
                    cell.begin(), cell.end(), [&](const Particle& p) {
                        if (domain.contains(p.coord))
                            return false;

                        for (const auto& cond : conditions_) {
                            cond->apply_to_particle(p, species_name, emitter,
                                                    domain);
                        }
                        return true;
                    });

                cell.erase(it, cell.end());
            }
        }
    }
    int count = 0;
    // добавить новые частицы текущего сорта
    for (const auto& p : emitter.current_species_particles()) {
        particles.add_particle(p);
        count++;
    }
    emitter.clear_current_species_buffer();
    std::cout << "Added " << count << " " << species_name
              << " particles from emitter\n";
    std::cout << "Time to apply boundary conditions: "
              << omp_get_wtime() - time0 << " seconds\n";
}

void BoundaryConditionHandler::flush_species(
    std::unordered_map<std::string, std::unique_ptr<ParticlesArray>>&
        all_species) {
    // добавить частицы в другие сорта по имени
    const auto& other = emitter.other_species_particles();

    for (const auto& kv : other) {
        int count = 0;
        const std::string& name = kv.first;
        if (all_species.count(name) == 0) {
            throw std::runtime_error("Species " + name +
                                     " not found when adding from emitter");
        }
        auto& target_species = all_species[name];
        for (const auto& p : kv.second) {
            target_species->add_particle(p);
            count++;
        }
        std::cout << "Added " << count << " " << name
                  << " particles from emitter\n";
    }
    emitter.clear_other_species_buffers();
}

void OpenBoundaryCondition::apply_to_operator(Operator& mat, const Domain& domain) const{
    const auto size = domain.size();
    mat.makeCompressed();
    // Получаем указатели на внутренние данные CSR
    double* values = mat.valuePtr();   // Массив значений
    const int* outerIndex =
        mat.outerIndexPtr();   // Массив индексов начала строк
    const int* innerIndex = mat.innerIndexPtr();   // Массив индексов столбцов

    auto set_values_zero = [](double* values, Face face, const Domain& domain,
                              int vindg, int vindg1, int value_index) {
        int i = domain.grid.pos_vind(vindg, 0);
        int j = domain.grid.pos_vind(vindg, 1);
        int k = domain.grid.pos_vind(vindg, 2);
        int d = domain.grid.pos_vind(vindg, 3);
        int i1 = domain.grid.pos_vind(vindg1, 0);
        int j1 = domain.grid.pos_vind(vindg1, 1);
        int k1 = domain.grid.pos_vind(vindg1, 2);
        int d1 = domain.grid.pos_vind(vindg1, 3);
        auto pos1 = domain.get_node_position(i, j, k, FieldType::ELECTRIC, d);
        auto pos2 = domain.get_node_position(i1, j1, k1, FieldType::ELECTRIC, d1);
        bool setZero1 = domain.geom.is_outside_face(face, pos1, 1.e-12);
        bool setZero2 = domain.geom.is_outside_face(face, pos2, 1.e-12);
        if (setZero1 || setZero2) {
            values[value_index] = 0.;
        }
    };

#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < 3 * (size.x() * size.y() * size.z()); i++) {
        // int ix = domain.grid.pos_vind(i, 0);
        // int iy = domain.grid.pos_vind(i, 1);
        // int iz = domain.grid.pos_vind(i, 2);
        // int d = domain.grid.pos_vind(i, 3);
        // auto pos = domain.get_node_position(ix, iy, iz, FieldType::ELECTRIC, d);
        // bool setZero = domain.geom.is_outside_face(face, pos, 1.e-12);
        // if (!domain.geom.is_outside_face(face_, pos, 1.e-12))
        //     continue;

        // Получаем диапазон ненулевых элементов для строки i
        int rowStart = outerIndex[i];
        int rowEnd = outerIndex[i + 1];

        for (int j = rowStart; j < rowEnd; j++) {
            int col = innerIndex[j];   // Столбец текущего элемента
            set_values_zero(values, face_, domain, i, col, j);
        }
    }

}
