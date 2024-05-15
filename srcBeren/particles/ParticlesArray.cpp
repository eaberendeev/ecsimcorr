#include "containers.h"
#include "World.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "bounds.h"
#include "service.h"
#include "collision.h"

//int ParticlesArray::counter;
ParticlesArray::ParticlesArray(const ParametersMap& particlesParams,
                               const ParametersMap& parameters,
                               const Domain& domain)
    : particlesData(domain.size()),
      countInCell(domain.size()),
      densityOnGrid(domain.size()),
      phaseOnGrid(100, 100), // TO DO: delete magick
      currentOnGrid(domain.size(), 3),
      Pxx(domain.size().x(), domain.size().y(), domain.size().z()),
      Pyy(domain.size().x(), domain.size().y(), domain.size().z()),
      Pzz(domain.size().x(), domain.size().y(), domain.size().z()),
      xCellSize(domain.cell_size().x()),
      yCellSize(domain.cell_size().y()),
      zCellSize(domain.cell_size().z()),
      xCellCount(domain.num_cells().x()),
      yCellCount(domain.num_cells().y()),
      zCellCount(domain.num_cells().z()),
      _name(particlesParams.get_values("Particles").at(0)),
      temperature(particlesParams.get_double("Temperature")),
      _mass(particlesParams.get_double("Mass")),
      charge(particlesParams.get_double("Charge")),
      density(particlesParams.get_double("Density")),
      NumPartPerCell(parameters.get_int("NumPartPerCell")),
      _mpw(density / NumPartPerCell) {
    injectionEnergy = lostEnergy = 0.;
}

// TO DO: change vector of particles to map of particles (key is name)
int get_num_of_type_particles(const std::vector<ParticlesArray>& Particles,
                              const std::string& ParticlesType) {
    for (int i = 0; i < 10; ++i) {
        if (Particles[i].name() == ParticlesType)
            return i;
    }
    return -1;
}

void ParticlesArray::add_particle(const Particle &particle){
    
    double xk = particle.coord.x() / xCellSize + GHOST_CELLS;
    double yk = particle.coord.y() / yCellSize + GHOST_CELLS;
    double zk = particle.coord.z() / zCellSize + GHOST_CELLS;

    int ix = int(xk);
    int iy = int(yk);
    int iz = int(zk);
    particlesData(ix,iy,iz).push_back(particle);
}

void ParticlesArray::save_init_coord_and_velocity() {
#pragma omp parallel for
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){
            particle.initCoord = particle.coord;
            particle.initVelocity = particle.velocity;
        }
    }
}
void ParticlesArray::save_init_coord() {
#pragma omp parallel for
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.initCoord = particle.coord;
        }
    }
}
void ParticlesArray::save_init_velocity() {
#pragma omp parallel for
    for (auto k = 0; k < size(); ++k) {
        for (auto& particle : particlesData(k)) {
            particle.initVelocity = particle.velocity;
        }
    }
}

// TO DO - move to simulation class
void ParticlesArray::delete_bounds(){
    lostEnergy = 0;
    Particle particle;
    for ( int ix = 0; ix < particlesData.size().x(); ++ix){
		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
		    	int ip = 0;
            	while( ip < countInCell(ix,iy,iz) ){  

                    if(ix > 8 && ix < particlesData.size().x() - 9 &&
                     iy > 8 && iy < particlesData.size().y() - 9 ){
                     //iz > 8 && iz < particlesData.size().z() - 9){
					    ip++;
				    } else {
					    delete_particle_runtime(ix,iy,iz,ip);
                                            lostEnergy += get_energy_particle(
                                                particle.velocity, _mass, _mpw);
                        }
                }
            }
        }
    }
    update_count_in_cell();
}

void ParticlesArray::update_cells(const Domain& domain){

    Particle particle;
    int ix2,iy2,iz2;
    for ( int ix = 0; ix < particlesData.size().x(); ++ix){
		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
		    	int ip = 0;
            	while( ip < countInCell(ix,iy,iz) ){  
                    particle = particlesData(ix,iy,iz)[ip];  
                    ix2 = int(particle.coord.x() / xCellSize + GHOST_CELLS);
                    iy2 = int(particle.coord.y() / yCellSize + GHOST_CELLS);
                    iz2 = int(particle.coord.z() / zCellSize + GHOST_CELLS);
                    if(ix == ix2 && iy == iy2 && iz == iz2){
                        particlesData(ix,iy,iz)[ip] = particle;  
					    ip++;
				    } else {
					    delete_particle_runtime(ix,iy,iz,ip);
                        if( particle_boundaries(particle.coord, domain) ){
                            ix2 = int(particle.coord.x() / xCellSize +
                                      GHOST_CELLS);
                            iy2 = int(particle.coord.y() / yCellSize +
                                      GHOST_CELLS);
                            iz2 = int(particle.coord.z() / zCellSize +
                                      GHOST_CELLS);
                            particlesData(ix2,iy2,iz2).push_back(particle);
                        } else{ 
                            lostEnergy += 0.5 * _mpw * _mass * dot(particle.velocity,particle.velocity); 
                        }
                    }
                }
            }
        }
    }
    update_count_in_cell();
}

void ParticlesArray::prepare(int timestep){
    currentOnGrid.set_zero();
    save_init_coord_and_velocity();
}

bool ParticlesArray::particle_boundaries(double3& coord, const Domain& domain)
{
    double sizeX = domain.cell_size(Dim::X)* domain.num_cells(Dim::X);
    double sizeY = domain.cell_size(Dim::Y) * domain.num_cells(Dim::Y);
    double sizeZ = domain.cell_size(Dim::Z) * domain.num_cells(Dim::Z);

    bool lostX = (coord.x() < 0. || coord.x() >= sizeX );
    bool lostY = (coord.y() < 0. || coord.y() >= sizeY );
    bool lostZ = (coord.z() < 0. || coord.z() >= sizeZ );
                    
    bool inArea = ! (lostX || lostY || lostZ);
      if (!inArea){
        if(isPeriodicX){
            if( coord.x() < 0. ){
                coord.x() += sizeX;
            }
            if( coord.x() >= sizeX ){
                coord.x() -= sizeX;
            }
        }
        if(isPeriodicY){
            if( coord.y() < 0. ){
                coord.y() += sizeY;
            }
            if( coord.y() >= sizeY ){
                coord.y() -= sizeY;
            }
        }
        if(isPeriodicZ){
            if( coord.z() < 0. ){
                coord.z() += sizeZ;
            }
            if( coord.z() >= sizeZ ){
                coord.z() -= sizeZ;
            }
        }
      }

                    
    lostX = (coord.x() < 0. || coord.x() >= sizeX );
    lostY = (coord.y() < 0. || coord.y() >= sizeY );
    lostZ = (coord.z() < 0. || coord.z() >= sizeZ );
                    
    return  !(lostX || lostY || lostZ);
}
