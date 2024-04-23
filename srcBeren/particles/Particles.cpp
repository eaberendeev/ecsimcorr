#include "containers.h"
#include "World.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "bounds.h"
#include "service.h"
#include "collision.h"
#include "recovery.h"

int ParticlesArray::counter;
ParticlesArray::ParticlesArray(const std::vector<std::string>& vecStringParams, World& world, const Domain &domain) :
    particlesData(world.region.numNodes),
    countInCell(world.region.numNodes),
    densityOnGrid(world.region.numNodes),phaseOnGrid(PxMax,PpMax),
    currentOnGrid(world.region.numNodes,3),
    Pxx(world.region.numNodes(0),world.region.numNodes(1),world.region.numNodes(2)),
    Pyy(world.region.numNodes(0),world.region.numNodes(1),world.region.numNodes(2)),
    Pzz(world.region.numNodes(0),world.region.numNodes(1),world.region.numNodes(2)),
    xCellSize(domain.cell_size().x()), yCellSize(domain.cell_size().y()), zCellSize(domain.cell_size().z()),
    xCellCount(domain.num_cells().x()), yCellCount(domain.num_cells().y()), zCellCount(domain.num_cells().z())
{
    for (const auto& line: vecStringParams){
        set_params_from_string(line);
    }
    injectionEnergy = lostEnergy =  0.;

}

int get_num_of_type_particles(const std::vector<ParticlesArray> &Particles, const std::string& ParticlesType){
  for (int i = 0; i < NumOfPartSpecies; ++i){
    if(Particles[i].name() == ParticlesType) return i;
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

// Устанавливаем значения основных параметров в переменной Params в соответствии с содержимым сроки
void ParticlesArray::set_params_from_string(const std::string& line){
    std::vector<std::string> strvec;

    strvec = split(line, ' ');

    if(strvec[0]=="Particles") {
        this->_name = strvec[1];
    }

    if(strvec[0]=="Temperature"){
        temperature = stod(strvec[1]);
    }

    if(strvec[0]=="Mass"){
        this->_mass = stod(strvec[1]);
    }

    if(strvec[0]=="Charge"){
       this->charge =  stod(strvec[1]);
    }
    if(strvec[0]=="SmoothMass"){
       option.smoothMass =  stol(strvec[1]);
    }
    if(strvec[0]=="SmoothMassSize"){
       option.smoothMassSize =  stod(strvec[1]);
    }
    if(strvec[0]=="SmoothMassMax"){
       option.smoothMassMax =  stod(strvec[1]);
    }  
    if(strvec[0]=="WidthY"){
        widthY = stod(strvec[1]);
    }
    if(strvec[0]=="WidthZ"){
        widthZ = stod(strvec[1]);
    }
    if(strvec[0]=="Density"){
        this->density = stod(strvec[1]); 
        this->_mpw = density / double(NumPartPerCell);
    }
    if(strvec[0]=="Focus"){
       focus = stod(strvec[1]);
    }
    if(strvec[0]=="Velocity"){
        velocity = stod(strvec[1]);
    }
    if(strvec[0]=="Pot_I"){
        this->pot_I = stod(strvec[1]);
    }
    if(strvec[0]=="Pot_k"){
        this->pot_k = stod(strvec[1]);
    }

    if(strvec[0]=="Px_min"){
        this->phasePXmin = stod(strvec[1]);
    }

    if(strvec[0]=="Px_max"){
         this->phasePXmax = stod(strvec[1]);
    }
    if (strvec[0]=="DistParams"){
        initDist = strvec[1];
        if(strvec[1]=="UniformCosX_dn_k"){
            distParams.push_back(stod(strvec[2]));
            distParams.push_back(stod(strvec[3]));
        }
    }
    if(strvec[0]=="BoundResumption"){
        option.boundResumption = stod(strvec[1]);
    }
    if(strvec[0]=="SourceType"){
        option.sourceType = stod(strvec[1]);
    }
    if(strvec[0]=="sourceAngle"){
        option.sourceAngle = stod(strvec[1]);
    }
}

void ParticlesArray::save_init_coord() {
#pragma omp parallel for
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){
            particle.initCoord = particle.coord;
            particle.initVelocity = particle.velocity;
        }
    }
}

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

    int counter = 0;
    for(auto cell = 0; cell < size(); ++cell){
        counter += particlesData(cell).size();
    }
    std::cout << "number of particles = " << counter << "\n";
    save_init_coord();
    
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
