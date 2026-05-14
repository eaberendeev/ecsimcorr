#include "ParticlesArray.h"

#include "Shape.h"
#include "World.h"
#include "collision.h"
#include "containers.h"

ParticlesArray::ParticlesArray(const nlohmann::json& config, const Domain& domain)
    : particlesData(domain.size()),
      densityOnGrid(domain.size(), 1),
      currentOnGrid(domain.size(), 3),
      charge(config.value("Charge", 1.0)),
      density(config.value("Density", 1.0)),
      name_(config.value("Name", "")),
      mass_(config.value("Mass", 1.0)),
      domain_(domain),
      config_(config) {
    NumPartPerCell = config_["NumPartPerCell"].get<int>();
    mpw_ = (density / NumPartPerCell);

    injectionEnergy = lostEnergyZ = lostEnergyXY = 0.;
    lostParticlesXY = lostParticlesZ = 0;

    initialize_distributions(config_);
}

ParticlesArray* find_species(Species& species, const std::string& name) {
    auto it = species.find(name);
    return (it == species.end()) ? nullptr : it->second.get();
}

const ParticlesArray* find_species(const Species& species, const std::string& name) {
    auto it = species.find(name);
    return (it == species.end()) ? nullptr : it->second.get();
}

void ParticlesArray::add_particle(const Particle& particle) {
    Vector3I cell_id = domain_.get_cell_index(particle.coord);

    particlesData(cell_id.x(), cell_id.y(), cell_id.z()).push_back(particle);
}

void ParticlesArray::add_particles(const std::vector<Particle>& particles) {
    for (const auto& particle : particles) {
        add_particle(particle);
    }
}

void ParticlesArray::save_init_coord_and_velocity() {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.initCoord = particle.coord;
            particle.initVelocity = particle.velocity;
        }
    }
}
void ParticlesArray::save_init_coord() {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.initCoord = particle.coord;
        }
    }
}
void ParticlesArray::save_init_velocity() {
#pragma omp parallel for schedule(dynamic, 32)
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.initVelocity = particle.velocity;
        }
    }
}

void ParticlesArray::update_cells(const Domain& domain) {
    const int COLOR_DIV = 3;
    const int COLOR_COUNT = 27;   // 3^3

    const int nx = particlesData.size().x();
    const int ny = particlesData.size().y();
    const int nz = particlesData.size().z();

    auto is_cell_of_color = [&](int ix, int iy, int iz, int color) {
        int cell_color = (ix % COLOR_DIV) + COLOR_DIV * ((iy % COLOR_DIV) + COLOR_DIV * (iz % COLOR_DIV));
        return cell_color == color;
    };

#pragma omp parallel
    {
        for (int color = 0; color < COLOR_COUNT; ++color) {
#pragma omp for schedule(static) collapse(3)
            for (int ix = 0; ix < nx; ++ix) {
                for (int iy = 0; iy < ny; ++iy) {
                    for (int iz = 0; iz < nz; ++iz) {
                        if (!is_cell_of_color(ix, iy, iz, color))
                            continue;

                        auto& cell_particles = particlesData(ix, iy, iz);
                        int ip = 0;
                        while (ip < static_cast<int>(cell_particles.size())) {
                            Particle particle = cell_particles[ip];
                            const Vector3I cell_id = domain.get_cell_index(particle.coord);

                            const int dx = std::abs(cell_id.x() - ix);
                            const int dy = std::abs(cell_id.y() - iy);
                            const int dz = std::abs(cell_id.z() - iz);
                            if (dx > 1 || dy > 1 || dz > 1) {
                                std::cout << "particle move error " << particle << "\n";
                                exit(0);
                            }

                            if (cell_id == Vector3I(ix, iy, iz)) {
                                ++ip;
                            } else {
                                std::swap(cell_particles[ip], cell_particles.back());
                                cell_particles.pop_back();
                                auto [ix2, iy2, iz2] = cell_id.split();
                                particlesData(ix2, iy2, iz2).push_back(particle);
                            }
                        }
                    }
                }
            }   // sync
        }
    }
}

// // Check periodic boundaries once and store result
// const bool has_periodic_bound = domain.is_periodic_bound(Axis::X) ||
//                                 domain.is_periodic_bound(Axis::Y) ||
//                                 domain.is_periodic_bound(Axis::Z);

// if (has_periodic_bound) {
//     for (int ix = 0; ix < nx; ++ix) {
//         for (int iy = 0; iy < ny; ++iy) {
//             for (int iz = 0; iz < nz; ++iz) {
//                 if (!domain.is_ghost_cell(ix, iy, iz))
//                     continue;

//                 for (const auto& original_particle :
//                      particlesData(ix, iy, iz)) {
//                     auto periodic_particle = original_particle;
//                     domain.make_point_periodic(periodic_particle.coord);
//                     add_particle(periodic_particle);
//                 }
//                 particlesData(ix, iy, iz).clear();
//             }
//         }
//     }
// }

// }

// void ParticlesArray::make_second_emisson(const Particle& particle){
//     auto [isInside, axis] = domain.in_region(particle.coord);
//     if (axis == Axis::Z){
//         add_particle(new_particle);
//     }
// }

void ParticlesArray::prepare() {
    currentOnGrid.setZero();
    save_init_coord_and_velocity();
}
