#include "containers.h"
#include "World.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "service.h"
#include "collision.h"

//int ParticlesArray::counter;
ParticlesArray::ParticlesArray(const ParametersMap& particlesParams,
                               const ParametersMap& parameters,
                               const Domain& domain)
    : particlesData(domain.size()),
      countInCell(domain.size()),
      densityOnGrid(domain.size(), 1),
      currentOnGrid(domain.size(), 3),
      charge(particlesParams.get_double("Charge")),
      density(particlesParams.get_double("Density")),
      _name(particlesParams.get_values("Particles").at(0)),
      temperature(particlesParams.get_double("Temperature", 0),
                  particlesParams.get_double("Temperature", 1),
                  particlesParams.get_double("Temperature", 2)),
      NumPartPerCell(parameters.get_int("NumPartPerCell")),

      _mass(particlesParams.get_double("Mass")),
      _mpw(density / NumPartPerCell),
      xCellSize(domain.cell_size().x()),
      yCellSize(domain.cell_size().y()),
      zCellSize(domain.cell_size().z()),
      xCellCount(domain.num_cells().x()),
      yCellCount(domain.num_cells().y()),
      zCellCount(domain.num_cells().z()),
      sortParameters(particlesParams) {
    countInCell.setZero();
    injectionEnergy = lostEnergyZ = lostEnergyXY = 0.;
    lostParticlesXY = lostParticlesZ = 0;

    distType = particlesParams.get_values("DistType").at(0);
    distSpace = particlesParams.get_values("DistSpace");
    distPulse = particlesParams.get_values("DistPulse");

#ifdef SET_PARTICLE_IDS
    ids = 0;
#endif
    // if(distType == "INITIAL"){
    //     distribute_particles(parameters,0);
    // }
    bounds = domain.get_bounds();
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

#ifdef SET_PARTICLE_IDS
    particle.id = ids++;
#endif
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

    // periodic boundary conditions
    if (domain.is_periodic_bound(Dim::X) || domain.is_periodic_bound(Dim::Y) || domain.is_periodic_bound(Dim::Z) )
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                if (!domain.is_ghost_cell(ix,iy,iz)) continue;
                for (auto& particle : particlesData(ix, iy, iz)) {
                    domain.make_point_periodic(particle.coord);
                    auto [ix2, iy2, iz2] =
                        get_cell_index(particle.coord);
                    particlesData(ix2, iy2, iz2).push_back(particle);
                }
                particlesData(ix, iy, iz).clear();
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

bool ParticlesArray::particle_boundaries(Particle& particle, const Domain& domain) {
    if(is_neutral()){
        auto isInside = domain.check_bbox_dim_bool(particle.coord, -1);
        if (isInside) {
            return true;
        } else {
            return false;
        }
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
