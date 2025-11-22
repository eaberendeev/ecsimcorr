#include "containers.h"
#include "World.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "service.h"
#include "collision.h"

ParticlesArray::ParticlesArray(
                               const nlohmann::json& config,
                               const ParametersMap& parameters,
                               const Domain& domain)
    : particlesData(domain.size()),
      countInCell(domain.size()),
      densityOnGrid(domain.size(), 1),
      currentOnGrid(domain.size(), 3),
      charge(config.value("Charge", 1.0)),
      density(config.value("Density", 1.0)),
      _name(config.value("Name", "")),
      _mass(config.value("Mass", 1.0)),
      xCellSize(domain.cell_size().x()),
      yCellSize(domain.cell_size().y()),
      zCellSize(domain.cell_size().z()),
      xCellCount(domain.num_cells().x()),
      yCellCount(domain.num_cells().y()),
      zCellCount(domain.num_cells().z()),
      config(config) {
    countInCell.setZero();
    NumPartPerCell = config["NumPartPerCell"].get<int>();
    _mpw = (density / NumPartPerCell) ;

    injectionEnergy = lostEnergyZ = lostEnergyXY = 0.;
    lostParticlesXY = lostParticlesZ = 0;

    bounds = domain.get_bounds();
    initializeDistributions(config);
}
// TO DO: change vector of particles to map of particles (key is name)
int get_num_of_type_particles(const Species& species,
                              const std::string& ParticlesType) {
    for (size_t i = 0; i < species.size(); ++i) {
        if (species[i]->name() == ParticlesType)
            return i;
    }
    return -1;
}

void ParticlesArray::add_particle(Particle &particle){
    
    double xk = particle.coord.x() / xCellSize + GHOST_CELLS;
    double yk = particle.coord.y() / yCellSize + GHOST_CELLS;
    double zk = particle.coord.z() / zCellSize + GHOST_CELLS;

    int ix = int(xk);
    int iy = int(yk);
    int iz = int(zk);

    particlesData(ix,iy,iz).push_back(particle);
}

void ParticlesArray::add_particles(std::vector<Particle> &particles){
    for(auto& particle : particles){
        add_particle(particle);
    }
    update_count_in_cell();
}

void ParticlesArray::save_init_coord_and_velocity() {
#pragma omp parallel for schedule(dynamic, 32)
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){
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
        int cell_color =
            (ix % COLOR_DIV) +
            COLOR_DIV * ((iy % COLOR_DIV) + COLOR_DIV * (iz % COLOR_DIV));
        return cell_color == color;
    };

    for (int color = 0; color < COLOR_COUNT; ++color) {
#pragma omp parallel for collapse(3) schedule(dynamic)
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    if (!is_cell_of_color(ix, iy, iz, color))
                        continue;

                    int ip = 0;
                    while (ip < countInCell(ix, iy, iz)) {
                        Particle particle = particlesData(ix, iy, iz)[ip];
                        auto [ix2, iy2, iz2] = get_cell_index(particle.coord);
                        if (ix == ix2 && iy == iy2 && iz == iz2) {
                            ip++;
                        } else {
                            delete_particle_runtime(ix, iy, iz, ip);
                            if (particle_boundaries(particle, domain)) {
                                particlesData(ix2, iy2, iz2)
                                    .push_back(particle);
                            }
                        }
                    }
                }
            }
        }
    }

    // Check periodic boundaries once and store result
    const bool has_periodic_bound = domain.is_periodic_bound(Dim::X) ||
                                    domain.is_periodic_bound(Dim::Y) ||
                                    domain.is_periodic_bound(Dim::Z);

    if (has_periodic_bound) {
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int iz = 0; iz < nz; ++iz) {
                    if (!domain.is_ghost_cell(ix, iy, iz))
                        continue;

                    for (const auto& original_particle :
                         particlesData(ix, iy, iz)) {
                        auto periodic_particle = original_particle;
                        domain.make_point_periodic(periodic_particle.coord);
                        add_particle(periodic_particle);
                    }
                    particlesData(ix, iy, iz).clear();
                }
            }
        }
    }

    update_count_in_cell();
}

void ParticlesArray::prepare(){
    currentOnGrid.setZero();
    save_init_coord_and_velocity();
    update_count_in_cell();
}

bool ParticlesArray::particle_boundaries(Particle& particle,
                                         const Domain& domain) {
    if (is_neutral()) {
        return domain.check_bbox_dim_bool(particle.coord, -1);
    }

    auto [isInside, axis] = domain.in_region(particle.coord);
    if (!isInside) {
        const double energy =
            get_energy_particle(particle.velocity, _mass, _mpw);
        if (axis == Axis::Z) {
            lostEnergyZ += energy;
            lostParticlesZ++;
        } else {
            lostEnergyXY += energy;
            lostParticlesXY++;
        }
        return false;
    }
    return true;
}
