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

void PeriodicParticlesCondition::apply_to_particle(const Particle& p,
                                                   BoundaryEmitter& emitter,
                                                   const Domain& domain) const {
    Particle p_new = p;
    domain.make_point_periodic(p_new.coord);
    emitter.emit_current_species(p_new);
}

void SecondEmissionCondition::apply_to_particle(const Particle& p,
                                                 BoundaryEmitter& emitter,
                                                 const Domain& domain) const {
//     // if(sort_name != "Electrons") return;
//     Particle p_new = p;
//     // p_new.coord[dim] = is_lower? -p_new.coord[dim] : p_new.coord[dim] - domain.cell_size()[dim]*domain.num_cells()[dim];

//     emitter.emit_to_species("Ions",p_new);
}

void BoundaryConditionHandler::apply_to_particles(
    ParticlesArray& particles,
    std::vector<std::unique_ptr<ParticlesArray>>& all_species,
    const Domain& domain) const {
    double time0 = omp_get_wtime();
    auto& data = particles.particlesData;
    int nx = data.size().x(), ny = data.size().y(), nz = data.size().z();

    BoundaryEmitter emitter(all_species.size());

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                // if (!domain.is_cell_outside(ix, iy, iz))
                //     continue;

                auto& cell = data(ix, iy, iz);

                auto it = std::remove_if(
                    cell.begin(), cell.end(), [&](const Particle& p) {
                        if (domain.in_region(p.coord, -1))
                            return false;

                        for (const auto& cond : conditions_) {
                            cond->apply_to_particle(p, emitter, domain);
                        }

                        return true;
                    });

                cell.erase(it, cell.end());
            }
        }
    }
    // добавить новые частицы текущего сорта

    for (const auto& p : emitter.current_species_particles()){
        particles.add_particle(p);
    }
    // добавить частицы в другие сорта

    const auto& other = emitter.other_species_particles();

    for (size_t s = 0; s < other.size(); ++s) {
        for (const auto& p : other[s]) {all_species[s]->add_particle(p);}
    }

    std::cout << "Time to apply boundary conditions: "
              << omp_get_wtime() - time0 << " seconds\n";
}
