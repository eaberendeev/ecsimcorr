#include "Vec.h"
#include "World.h"
#include "ParticlesArray.h"
#include "Shape.h"
#include "bounds.h"
#include "service.h"
#include "collision.h"

std::ostream& operator<<(std::ostream& out,const double3& val){
	  out <<  val.x() << " " << val.y() << " " << val.z();
	  return out;
} 

std::ostream& operator<<(std::ostream& out, const ParticleSimple& particle){
	out << particle.coord << " " << particle.velocity;
	return out;
} 
std::ostream& operator<<(std::ostream& out, const ParticleMass& particle){
	out << particle.coord << " " << particle.velocity << " " << particle.mass;
	return out;
} 




int ParticlesArray::counter;
ParticlesArray::ParticlesArray(const std::vector<std::string>& vecStringParams, World& world):
    particlesData(world.region.numNodes),
    countInCell(world.region.numNodes),
    densityOnGrid(world.region.numNodes),phaseOnGrid(PxMax,PpMax),
    currentOnGrid(world.region.numNodes,3),
    Pxx(world.region.numNodes(0),world.region.numNodes(1),world.region.numNodes(2)),
    Pyy(world.region.numNodes(0),world.region.numNodes(1),world.region.numNodes(2)),
    Pzz(world.region.numNodes(0),world.region.numNodes(1),world.region.numNodes(2)),
    _world(world) 
{
    for (const auto& line: vecStringParams){
        set_params_from_string(line);
    }
    injectionEnergy = lostEnergy =  0.;

    if (RECOVERY > 0) {
        read_particles_from_recovery();
        std::cout << "Upload " + name + " success!\n";
    }
    else {
        if( k_particles_reservation > 0. ){
            for(auto k = 0; k < size(); ++k){
                    particlesData(k).reserve(int(k_particles_reservation*NumPartPerCell));
            }
        }
        set_distribution();
    }
    update_count_in_cell();
    density_on_grid_update();
    std::cout << particlesData.size() << " " << particlesData.capacity()  <<"\n";

}

double PulseFromKev(double kev, double mass){
  double gama = kev / MC2 + mass;
  return sqrt((gama*gama)- mass);
}


int get_num_of_type_particles(const std::vector<ParticlesArray> &Particles, const std::string& ParticlesType){
  for (int i = 0; i < NumOfPartSpecies; ++i){
    if(Particles[i].name == ParticlesType) return i;
  }
  return -1;
}

void ParticlesArray::write_particles_to_recovery(int timestep){
    if( RecTimeStep < 0){
        return;
    }
    if( timestep % RecTimeStep != 0){
        return;
    }

    std::cout << "Backup " + name + " in " << timestep << "\n";
    Particle particle;
    std::ofstream file_bin(".//Recovery//Particles//"+name+"//"+name+".backup", std::ios::out | std::ios::binary);
    int n_particles = 0;
    for(auto cell = 0; cell < size(); ++cell){
            n_particles += particlesData(cell).size();   
    }

    file_bin.write((char*) &n_particles, sizeof(n_particles));
    for(auto cell = 0; cell < size(); ++cell){
        int cellSize = particlesData(cell).size();
        if(cellSize > 0){
            file_bin.write((char*) &particlesData(cell)[0], cellSize*sizeof(particlesData(cell)[0]));
        }
    }
    file_bin.close();
}

void ParticlesArray::read_particles_from_recovery(){

    Particle particle;
    std::ifstream file_bin("..//Recovery//Particles//" + name+"//"+ name+".backup", std::ios::in | std::ios::binary);
    int n_particles;
    file_bin.read((char*) &n_particles, sizeof(n_particles));
    for(int i = 0; i < n_particles; i++){
        file_bin.read((char*) &particle, sizeof(Particle));
        int ix = int(particle.coord.x() / Dx + CELLS_SHIFT);
        int iy = int(particle.coord.y() / Dy + CELLS_SHIFT);
        int iz = int(particle.coord.z() / Dz + CELLS_SHIFT);
        particlesData(ix,iy,iz).push_back(particle);
    }
    file_bin.close();
}

void ParticlesArray::add_particle(const Particle &particle){
    
    double xk = particle.coord.x() / Dx + CELLS_SHIFT;
    double yk = particle.coord.y() / Dy + CELLS_SHIFT;
    double zk = particle.coord.z() / Dz + CELLS_SHIFT;
    
    int ix = int(xk);
    int iy = int(yk);
    int iz = int(zk);
    particlesData(ix,iy,iz).push_back(particle);
}


void ParticlesArray::glue_density_bound()
{ 
    
    int overlap = 3;
    int3 size = densityOnGrid.size();
    
    if(isPeriodicX){
        for (auto i = 0; i < overlap; ++i){
          for (auto l = 0; l < size.y(); ++l){
            for (auto k = 0; k < size.z(); ++k){
                densityOnGrid(i,l,k) += densityOnGrid(i + size.x() - overlap,l,k);
                densityOnGrid(i + size.x() - overlap,l,k) = densityOnGrid(i,l,k);
            }
          }
        }
    }
    if(isPeriodicY){
        for (auto i = 0; i < size.x(); ++i){
          for (auto l = 0; l < overlap; ++l){
            for (auto k = 0; k < size.z(); ++k){
                densityOnGrid(i,l,k) += densityOnGrid(i,l + size.y() - overlap,k);
                densityOnGrid(i,l+ size.y() - overlap,k) = densityOnGrid(i,l ,k);
            }
          }
        }
    }
    if(isPeriodicZ){
        for (auto i = 0; i < size.x(); ++i){
          for (auto l = 0; l < size.y(); ++l){
            for (auto k = 0; k < overlap; ++k){
                densityOnGrid(i,l,k) += densityOnGrid(i ,l,k+ size.z() - overlap);
                densityOnGrid(i ,l,k+ size.z() - overlap) = densityOnGrid(i,l,k);
            }
          }
        }
    }
}

void ParticlesArray::density_on_grid_update(){ 
    constexpr auto SMAX = 2*SHAPE_SIZE;

    densityOnGrid.clear();

#pragma omp parallel for
    for ( auto j = 0; j < size(); ++j){
        double arg;
    alignas(64) double sx[SMAX];
    alignas(64) double sy[SMAX];
    alignas(64) double sz[SMAX];

        for(const auto& particle : particlesData(j)){
            int xk = int(particle.coord.x() / Dx);
            int yk = int(particle.coord.y() / Dy);
            int zk = int(particle.coord.z() / Dz);

            for(int n = 0; n < SMAX; ++n){
                arg = -particle.coord.y() / Dy + double(yk - CELLS_SHIFT + n);
                sy[n] = Shape2(arg) ;
                arg = -particle.coord.x() / Dx + double(xk - CELLS_SHIFT + n);
                sx[n] = Shape2(arg);
                arg = -particle.coord.z() / Dz + double(zk - CELLS_SHIFT + n);
                sz[n] = Shape2(arg);
            }

            for(int n = 0; n < SMAX; ++n){
                int indx = xk + n;
                for(int m = 0; m < SMAX; ++m){
                    int indy = yk + m;
                    for(int k = 0; k < SMAX; ++k){
                    int indz = zk + k; 
                    #pragma omp atomic update
                    densityOnGrid(indx,indy,indz) += _mpw*charge * sx[n] * sy[m] * sz[k];
                    } 
                }
            }
        }
    }
    glue_density_bound();

}

void ParticlesArray::phase_on_grid_update(){ 
    int pk, xk;
    double x, px;
    bool blounder;
    Particle particle;
    double x_min = Dx*(_world.region.dampCells[0].x());
    double x_max = Dx*(_world.region.numCells.x() - _world.region.dampCells[1].x() );
    double pdx = (x_max - x_min) / PxMax;
    double pdp = (phasePXmax -phasePXmin) / PpMax;

    phaseOnGrid.clear();

    for (auto k = 0 ; k < size(); k++){
        for(const auto& particle : particlesData(k)){
            x = particle.coord.x();
            px = particle.velocity.x();

            xk = int((x-x_min) / pdx); 
            pk = int((px - phasePXmin) / pdp);

            blounder = (xk<0) || (xk>=PxMax) || (pk<0) || (pk>=PpMax);
            if(!blounder){
                phaseOnGrid(xk,pk) += (mpw() / (Dx*Dy*pdx*pdp) );
            }
        }
    }
}


void ParticlesArray::set_distribution(){
		set_space_distribution();
    static ThreadRandomGenerator randomGenerator;
    set_pulse_distribution(randomGenerator);
}

// Устанавливаем значения основных параметров в переменной Params в соответствии с содержимым сроки
void ParticlesArray::set_params_from_string(const std::string& line){
    std::vector<std::string> strvec;

    strvec = split(line, ' ');

    if(strvec[0]=="Particles") {
        this->name = strvec[1];
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

double ParticlesArray::get_kinetic_energy() const{
    double energy = 0;
#pragma omp parallel for reduction(+:energy)
    for(auto k = 0; k < size(); ++k){
        for(const auto& particle : particlesData(k)){
            energy += get_energy_particle(particle.velocity, _mass, _mpw);
        }
    }

    return energy;
}

double ParticlesArray::get_kinetic_energy(int dim) const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            double3 velocity(0,0,0);
            velocity(dim) = particle.velocity(dim);
            energy +=
                get_energy_particle(velocity, _mass, _mpw);
        }
    }

    return energy;
}

double3 ParticlesArray::get_kinetic_energy_component() const {
    double enx = 0;
    double eny = 0;
    double enz = 0;
#pragma omp parallel for reduction(+ : enx) reduction(+ : eny) reduction(+ : enz)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            auto velocity = particle.velocity;
            enx += get_energy_particle(velocity.x(), _mass, _mpw);
            eny += get_energy_particle(velocity.y(), _mass, _mpw);
            enz += get_energy_particle(velocity.z(), _mass, _mpw);
        }
    }

    return double3(enx, eny, enz);
}
double ParticlesArray::get_kinetic_energy(int dim1, int dim2) const {
    double energy = 0;
#pragma omp parallel for reduction(+ : energy)
    for (auto k = 0; k < size(); ++k) {
        for (const auto& particle : particlesData(k)) {
            double3 velocity(0, 0, 0);
            velocity(dim1) = particle.velocity(dim1);
            velocity(dim2) = particle.velocity(dim2);
            energy += get_energy_particle(velocity, _mass, _mpw);
        }
    }

    return energy;
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

void ParticlesArray::update_cells(){

    Particle particle;
    int ix2,iy2,iz2;
    for ( int ix = 0; ix < particlesData.size().x(); ++ix){
		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
		    	int ip = 0;
            	while( ip < countInCell(ix,iy,iz) ){  
                    particle = particlesData(ix,iy,iz)[ip];  
                    ix2 = int(particle.coord.x() / Dx + CELLS_SHIFT);
                    iy2 = int(particle.coord.y() / Dy + CELLS_SHIFT);
                    iz2 = int(particle.coord.z() / Dz + CELLS_SHIFT);
                    if(ix == ix2 && iy == iy2 && iz == iz2){
                        particlesData(ix,iy,iz)[ip] = particle;  
					    ip++;
				    } else {
					    delete_particle_runtime(ix,iy,iz,ip);
                        if( particle_boundaries(particle.coord) ){
                            ix2 = int(particle.coord.x() / Dx + CELLS_SHIFT);
                            iy2 = int(particle.coord.y() / Dy + CELLS_SHIFT);
                            iz2 = int(particle.coord.z() / Dz + CELLS_SHIFT);
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

void ParticlesArray::move(double dt){
#pragma omp parallel for
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){        
                    particle.move(dt);
        }
    }
}

void ParticlesArray::prepare(int timestep){
    currentOnGrid.clear();

    int counter = 0;
    for(auto cell = 0; cell < size(); ++cell){
        counter += particlesData(cell).size();
    }
    std::cout << "number of particles = " << counter << "\n";
    save_init_coord();
    
}

bool ParticlesArray::particle_boundaries(double3& coord)
{
    double sizeX = Dx*_world.region.numCells.x();
    double sizeY = Dx*_world.region.numCells.y();
    double sizeZ = Dx*_world.region.numCells.z();
                                        
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

void ParticlesArray::predict_velocity(const Mesh& mesh){

#pragma omp parallel for
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){        
            const auto coord = particle.coord;
            const auto velocity = particle.velocity;
            double3 E, En, B;
            get_fields_in_pos(mesh.fieldE, mesh.fieldB,coord,E,B);
            En = get_fieldE_in_pos(mesh.fieldEp,coord); 
            E=0.5*(E+En);
            const auto beta = Dt * charge / _mass;
            const auto alpha = 0.5*beta*mag(B);
            const auto alpha2 = alpha*alpha;
            const auto h = unit(B);
            
            const auto v12 = (1./(1.+ alpha2 )) *
                                    (velocity + alpha*cross(velocity,h) + alpha2*dot(h,velocity)*h
                                    + 0.5*beta*( E + alpha*cross(E,h) + alpha2*dot(E,h)*h) );

            particle.velocity = 2.*v12 - velocity; 
        }   
    }
}

void ParticlesArray::move_and_calc_current(const double dt){

    constexpr auto SMAX = 2*SHAPE_SIZE;

    const double conx = Dx / (6*dt) * _mpw;
    const double cony = Dy / (6*dt) * _mpw;
    const double conz = Dz / (6*dt) * _mpw;

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        for(auto& particle : particlesData(pk)){        

                    double3 start = particle.coord;
                    
                    particle.move(dt);
                    
                    double3 end = particle.coord;

                    double xx = start.x() / Dx;
                    double yy = start.y() / Dy;
                    double zz = start.z() / Dz;

                    double xn = end.x() / Dx;
                    double yn = end.y() / Dy;
                    double zn = end.z() / Dz;
                    
                    int xk = int(xx);
                    int yk = int(yy);
                    int zk = int(zz);

                    for(int n = 0; n < SMAX; ++n){
                        for(int m = 0; m < SMAX; ++m){
                            for(int k = 0; k < SMAX; ++k){
                                jx[n][m][k] = 0.;
                                jy[n][m][k] = 0.;
                                jz[n][m][k] = 0.;
                            }
                        }
                    }
                    double arg;
                    for(int n = 0; n < SMAX; ++n){
                        arg = -xx + double(xk - CELLS_SHIFT + n);
                        sx[n] = Shape2(arg);
                        arg = -yy + double(yk - CELLS_SHIFT + n);
                        sy[n] = Shape2(arg);
                        arg = -zz + double(zk - CELLS_SHIFT + n);
                        sz[n] = Shape2(arg);            
                        arg = -xn + double(xk - CELLS_SHIFT + n);
                        sx_n[n] = Shape2(arg);
                        arg = -yn + double(yk - CELLS_SHIFT + n);
                        sy_n[n] = Shape2(arg);
                        arg = -zn + double(zk - CELLS_SHIFT + n);
                        sz_n[n] = Shape2(arg);
                    }

                    for(int n = 0; n < SMAX; ++n){
                        int indx = xk  + n;
                        for(int m = 0; m < SMAX; ++m){
                            int indy = yk + m;
                            for(int k = 0; k < SMAX; ++k){
                                int indz = zk + k;
                        
                                if(n == 0) jx[n][m][k] = -charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                                if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
                                
                                if(m == 0) jy[n][m][k] = -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                
                                if(k == 0) jz[n][m][k] = -charge * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-charge * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                            
                                #pragma omp atomic update
                                currentOnGrid(indx,indy,indz,0) += jx[n][m][k];
                                #pragma omp atomic update
                                currentOnGrid(indx,indy,indz,1) += jy[n][m][k];
                                #pragma omp atomic update
                                currentOnGrid(indx,indy,indz,2) += jz[n][m][k];
                            }
                        }
                    }
                }
            }

}

void ParticlesArray::get_P(){
    Pxx.clear();
    Pyy.clear();
    Pzz.clear();

    constexpr auto SMAX = SHAPE_SIZE;

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double wx[2];
        double wy[2];
        double wz[2];
        for(auto& particle : particlesData(pk)){        
                    const auto coord = particle.coord;
                    const auto velocity = particle.velocity;

                    auto x = coord.x() / Dx + CELLS_SHIFT;
                    auto y = coord.y() / Dy + CELLS_SHIFT;
                    auto z = coord.z() / Dz + CELLS_SHIFT;

                    const auto intx = int(x);
                    const auto inty = int(y);
                    const auto intz = int(z);

                    wx[1] = (x - intx);
                    wx[0] = 1 - wx[1];
                    wy[1] = (y - inty);
                    wy[0] = 1 - wy[1];
                    wz[1] = (z - intz);
                    wz[0] = 1 - wz[1];
                    
                    double vx = velocity.x();
                    double vy = velocity.y();
                    double vz = velocity.z();
                    double pxx = _mass*vx*vx * _mpw;
                    double pyy = _mass*vy*vy * _mpw;
                    double pzz = _mass*vz*vz * _mpw;


                    for(int nx = 0; nx < SMAX; ++nx){
                        const int i = intx + nx;
                        for(int ny = 0; ny < SMAX; ++ny){
                            const int j = inty  + ny;
                            for(int nz = 0; nz < SMAX; ++nz){
                                const int k = intz  + nz;
                                const auto sx = wx[nx] * wy[ny] * wz[nz];
                                #pragma omp atomic update
                                Pxx(i,j,k) += sx*pxx;          
                                #pragma omp atomic update
                                Pyy(i,j,k) += sx*pyy;          
                                #pragma omp atomic update
                                Pzz(i,j,k) += sx*pzz;
                            }
                        }
                    }

                }
            }
}

void ParticlesArray::get_Pr(){
    Pxx.clear();
    Pyy.clear();
    Pzz.clear();

    constexpr auto SMAX = SHAPE_SIZE;

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double wx[2];
        double wy[2];
        double wz[2];
        for (auto& particle : particlesData(pk)){

            const auto coord = particle.coord;
            const auto velocity = particle.velocity;

            auto x = coord.x() / Dx + CELLS_SHIFT;
            auto y = coord.y() / Dy + CELLS_SHIFT;
            auto z = coord.z() / Dz + CELLS_SHIFT;

            const auto intx = int(x);
            const auto inty = int(y);
            const auto intz = int(z);

            wx[1] = (x - intx);
            wx[0] = 1 - wx[1];
            wy[1] = (y - inty);
            wy[0] = 1 - wy[1];
            wz[1] = (z - intz);
            wz[0] = 1 - wz[1];

            double x0 = 0.5 * Dx * NumCellsX_glob;
            double R = sqrt((coord.x()-x0)*(coord.x()-x0) + (coord.y()-x0)*(coord.y()-x0) );

            double vr =   ((coord.x()-x0) / R) * velocity.x() + ((coord.y()-x0) / R) * velocity.y();
            double vp = - ((coord.y()-x0) / R) * velocity.x() + ((coord.x()-x0) / R) * velocity.y();
            double vz = velocity.z();
            double prr = _mass*vr*vr * _mpw;
            double ppp = _mass*vp*vp * _mpw;
            double pzz = _mass*vz*vz * _mpw;


            for(int nx = 0; nx < SMAX; ++nx){
                const int i = intx + nx;
               for(int ny = 0; ny < SMAX; ++ny){
                    const int j = inty  + ny;
                    for(int nz = 0; nz < SMAX; ++nz){
                        const int k = intz  + nz;
                        const auto sx = wx[nx] * wy[ny] * wz[nz];
                        #pragma omp atomic update
                        Pxx(i,j,k) += sx*prr; 
                        #pragma omp atomic update
                        Pyy(i,j,k) += sx*ppp;
                        #pragma omp atomic update
                        Pzz(i,j,k) += sx*pzz;     
                    } 
                }
            }
        }
    }

}

void ParticlesArray::move_and_calc_current(const double dt, Field3d& fieldJ){

    constexpr auto SMAX = 2*SHAPE_SIZE;

    const double conx = Dx / (6*dt) * _mpw;
    const double cony = Dy / (6*dt) * _mpw;
    const double conz = Dz / (6*dt) * _mpw;

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double arg;
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        for(auto& particle : particlesData(pk)){      
                    double3 start = particle.coord;
                    
                    particle.move(dt);
                    
                    double3 end = particle.coord;

                    double xx = start.x() / Dx;
                    double yy = start.y() / Dy;
                    double zz = start.z() / Dz;

                    double xn = end.x() / Dx;
                    double yn = end.y() / Dy;
                    double zn = end.z() / Dz;
                    
                    int xk = int(xx);
                    int yk = int(yy);
                    int zk = int(zz);

                    for(int n = 0; n < SMAX; ++n){
                        for(int m = 0; m < SMAX; ++m){
                            for(int k = 0; k < SMAX; ++k){
                                jx[n][m][k] = 0.;
                                jy[n][m][k] = 0.;
                                jz[n][m][k] = 0.;
                            }
                        }
                    }

                    for(int n = 0; n < SMAX; ++n){
                        arg = -xx + double(xk - CELLS_SHIFT + n);
                        sx[n] = Shape2(arg);
                        arg = -yy + double(yk - CELLS_SHIFT + n);
                        sy[n] = Shape2(arg);
                        arg = -zz + double(zk - CELLS_SHIFT + n);
                        sz[n] = Shape2(arg);            
                        arg = -xn + double(xk - CELLS_SHIFT + n);
                        sx_n[n] = Shape2(arg);
                        arg = -yn + double(yk - CELLS_SHIFT + n);
                        sy_n[n] = Shape2(arg);
                        arg = -zn + double(zk - CELLS_SHIFT + n);
                        sz_n[n] = Shape2(arg);
                    }

                    for(int n = 0; n < SMAX; ++n){
                        int indx = xk  + n;
                        for(int m = 0; m < SMAX; ++m){
                            int indy = yk + m;
                            for(int k = 0; k < SMAX; ++k){
                                int indz = zk + k;
                        
                                if(n == 0) jx[n][m][k] = -charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                                if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - charge * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
                                
                                if(m == 0) jy[n][m][k] = -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -charge * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                
                                if(k == 0) jz[n][m][k] = -charge * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-charge * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                            
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,0) += jx[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,1) += jy[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,2) += jz[n][m][k];
                            }
                        }
                    }
                }
            }
}

void ParticlesArray::calc_Esirkepov_current(const double dt, Field3d& fieldJ) const{
    constexpr auto SMAX = 2*SHAPE_SIZE;
        const double conx = Dx / (6*dt) * _mpw;
        const double cony = Dy / (6*dt) * _mpw;
        const double conz = Dz / (6*dt) * _mpw;
#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double arg;
        alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
        alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
        alignas(64) double jx[SMAX][SMAX][SMAX];
        alignas(64) double jy[SMAX][SMAX][SMAX];
        alignas(64) double jz[SMAX][SMAX][SMAX];

        double q = charge;

        for(auto& particle : particlesData(pk)){      
                    const auto start = particle.initCoord;                    
                    const auto end = particle.coord;

                    const auto xx = start.x() / Dx;
                    const auto yy = start.y() / Dy;
                    const auto zz = start.z() / Dz;

                    const auto xn = end.x() / Dx;
                    const auto yn = end.y() / Dy;
                    const auto zn = end.z() / Dz;
                    
                    const auto xk = int(xx);
                    const auto yk = int(yy);
                    const auto zk = int(zz);

                    for(int n = 0; n < SMAX; ++n){
                        for(int m = 0; m < SMAX; ++m){
                            for(int k = 0; k < SMAX; ++k){
                                jx[n][m][k] = 0.;
                                jy[n][m][k] = 0.;
                                jz[n][m][k] = 0.;
                            }
                        }
                    }

                    for(int n = 0; n < SMAX; ++n){
                        arg = -xx + double(xk - CELLS_SHIFT + n);
                        sx[n] = Shape(arg);
                        arg = -yy + double(yk - CELLS_SHIFT + n);
                        sy[n] = Shape(arg);
                        arg = -zz + double(zk - CELLS_SHIFT + n);
                        sz[n] = Shape(arg);            
                        arg = -xn + double(xk - CELLS_SHIFT + n);
                        sx_n[n] = Shape(arg);
                        arg = -yn + double(yk - CELLS_SHIFT + n);
                        sy_n[n] = Shape(arg);
                        arg = -zn + double(zk - CELLS_SHIFT + n);
                        sz_n[n] = Shape(arg);
                    }

                    for(int n = 0; n < SMAX; ++n){
                        const auto indx = xk  + n;
                        for(int m = 0; m < SMAX; ++m){
                            const auto indy = yk + m;
                            for(int k = 0; k < SMAX; ++k){
                                const auto indz = zk + k;
                        
                                if(n == 0) jx[n][m][k] = -q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

                                if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
                                
                                if(m == 0) jy[n][m][k] = -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                
                                if(k == 0) jz[n][m][k] = -q * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-q * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
                                #pragma omp atomic update                            
                                fieldJ(indx,indy,indz,0) += jx[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,1) += jy[n][m][k];
                                #pragma omp atomic update
                                fieldJ(indx,indy,indz,2) += jz[n][m][k];
                            }
                        }
                    }
                }
            }
}

void ParticlesArray::correctv(Mesh& mesh){

    std::array<double,20> ldistr;
    for (auto& val: ldistr){
        val =0.0;
    }

    double jp_cell = 0;
#pragma omp parallel for reduction(+:jp_cell)
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto initVelocity = particle.initVelocity;
                    const auto velocity = particle.velocity;
 
                    double3 end = particle.coord;
                    double3 coord = end - 0.5*Dt*velocity;
                    
                    double3 Ep = get_fieldE_in_pos(mesh.fieldEp,coord); 
                    double3 E = get_fieldE_in_pos(mesh.fieldE,coord);

                    double3 v12 = 0.5*(velocity + initVelocity); 

                    jp_cell += 0.5*_mpw*charge*dot(v12,(Ep+E));
                }
            }

    mesh.fieldJp_full.data() = mesh.fieldJp.data() + mesh.Lmat2*(mesh.fieldE.data() + mesh.fieldEp.data())/Dt;

    const double energyJeEn = mesh.calc_JE(mesh.fieldEn,currentOnGrid);
    const double energyJeE = mesh.calc_JE(mesh.fieldE,currentOnGrid);
    const double energyJpEp = mesh.calc_JE(mesh.fieldEp,mesh.fieldJp_full);
    const double energyJpE = mesh.calc_JE(mesh.fieldE,mesh.fieldJp_full);
    const double energyK = get_kinetic_energy();
    const double lambda = sqrt(1 + Dt*( 0.5*(energyJeEn + energyJeE) - jp_cell) / energyK);

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto velocity = particle.velocity;

                    particle.velocity = lambda*velocity;
                }
            }

    const double energyK2 = get_kinetic_energy();
    std::cout << "lambda "<< lambda  << " " << lambda*lambda << " "<< energyK2-energyK << " " << 0.5*Dt*(energyJeEn + energyJeE - energyJpEp - energyJpE) << "\n";

}


void ParticlesArray::correctv_component(Mesh& mesh){
    double jp_cellx = 0;
    double jp_celly = 0;
    double jp_cellz = 0;
#pragma omp parallel for reduction(+ : jp_cellx) reduction(+ : jp_celly) reduction(+:jp_cellz)
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto initVelocity = particle.initVelocity;
                    const auto velocity = particle.velocity;
 
                    double3 end = particle.coord;
                    double3 coord = end - 0.5*Dt*velocity;
                    
                    double3 Ep = get_fieldE_in_pos(mesh.fieldEp,coord); 
                    double3 E = get_fieldE_in_pos(mesh.fieldE,coord);
                    E+= Ep;

                    double3 v12 = 0.5*(velocity + initVelocity);
                    double3 vE = double3(v12.x() * E.x(), v12.y() * E.y(),
                                         v12.z() * E.z() );

                    jp_cellx += (0.5 * _mpw * charge) * vE.x();
                    jp_celly += (0.5 * _mpw * charge) * vE.y();
                    jp_cellz += (0.5 * _mpw * charge) * vE.z();
        }
            }


    const double3 energyJeEn = mesh.calc_JE_component(mesh.fieldEn,currentOnGrid);
    const double3 energyJeE =
        mesh.calc_JE_component(mesh.fieldE, currentOnGrid);
    const double3 energyK = get_kinetic_energy_component();
    double3 lambda;
    lambda.x() = 
            sqrt(1 + Dt * (0.5 * (energyJeEn.x() + energyJeE.x()) - jp_cellx) /
                     energyK.x());
    lambda.y() =
          sqrt(1 + Dt * (0.5 * (energyJeEn.y() + energyJeE.y()) - jp_celly) /
                     energyK.y());
    lambda.z() =
         sqrt(1 + Dt * (0.5 * (energyJeEn.z() + energyJeE.z()) - jp_cellz) /
                     energyK.z());
    // double lambda2 =
    //     sqrt(1 + Dt *
    //                  (0.5 * (energyJeEn.x() + energyJeE.x() + energyJeEn.y() +
    //                          energyJeE.y() + energyJeEn.z() + energyJeE.z()) -
    //                   jp_cellx - jp_celly - jp_cellz) /
    //                  (energyK.x() + energyK.y() + energyK.z()));

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto velocity = particle.velocity;

                    particle.velocity.x() = lambda.x() * velocity.x();
                    particle.velocity.y() = lambda.y() * velocity.y();
                    particle.velocity.z() = lambda.z() * velocity.z();
                    //particle.velocity = lambda2 * velocity;
        }
            }

    //const double energyK2 = get_kinetic_energy();
    std::cout << "lambda "<< lambda << "\n";

}

void ParticlesArray::predict_current(const Field3d& fieldB, Field3d& fieldJ){
    constexpr auto SMAX = SHAPE_SIZE; 
#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        int i,j,k;
        int i05,j05,k05;
        double qp = charge;

        alignas(64) double wx[SMAX], wy[SMAX], wz[SMAX];
        alignas(64) double wx05[SMAX], wy05[SMAX], wz05[SMAX];

        for(auto& particle : particlesData(pk)){     
                    double3 coord = particle.coord;
                    double3 velocity = particle.velocity;

                    double x = coord.x() / Dx + CELLS_SHIFT;
                    double y = coord.y() / Dy + CELLS_SHIFT;
                    double z = coord.z() / Dz + CELLS_SHIFT;
                    double x05 = x - 0.5;
                    double y05 = y - 0.5;
                    double z05 = z - 0.5;

                    const auto ix = int(x);
                    const auto iy = int(y);
                    const auto iz = int(z);
                    const auto ix05 = int(x05);
                    const auto iy05 = int(y05);
                    const auto iz05 = int(z05);
                    
                    wx[1] = (x - ix);
                    wx[0] = 1 - wx[1];
                    wy[1] = (y - iy);
                    wy[0] = 1 - wy[1];
                    wz[1] = (z - iz);
                    wz[0] = 1 - wz[1];

                    wx05[1] = (x05 - ix05);
                    wx05[0] = 1 - wx05[1];
                    wy05[1] = (y05 - iy05);
                    wy05[0] = 1 - wy05[1];
                    wz05[1] = (z05 - iz05);
                    wz05[0] = 1 - wz05[1];        
                    
                    double3 B = get_fieldB_in_pos(fieldB,coord); 

                    double beta = Dt * qp / _mass;
                    double alpha = 0.5*beta*mag(B);
                    double alpha2 = alpha*alpha;
                    double3 h = unit(B);

                    double3 current = qp * _mpw / (1.+alpha2) * (velocity + alpha*cross(velocity,h) + alpha2*dot(h,velocity)*h );

                    for(int nx = 0; nx < SMAX; ++nx){
                        i = ix + nx;
                        i05 = ix05 + nx;
                        for(int ny = 0; ny < SMAX; ++ny){
                            j = iy  + ny;
                            j05 = iy05  + ny;
                            for(int nz = 0; nz < SMAX; ++nz){
                                k = iz  + nz;
                                k05 = iz05  + nz;
                                double sx = wx05[nx] * wy[ny] * wz[nz];
                                double sy = wx[nx] * wy05[ny] * wz[nz];
                                double sz = wx[nx] * wy[ny] * wz05[nz];
                                #pragma omp atomic update
                                fieldJ(i05,j,k,0) += sx*current.x();
                                #pragma omp atomic update
                                fieldJ(i,j05,k,1) += sy*current.y();
                                #pragma omp atomic update
                                fieldJ(i,j,k05,2) += sz*current.z();                   
                            }
                        }
                    }

                }
            }
}

void ParticlesArray::get_L(Mesh& mesh) {
#pragma omp parallel
{
    for (int xStep = 0; xStep < 4; xStep++) {
        for (int yStep = 0; yStep < 4; yStep++) {
            #pragma omp for collapse(2)
            for (int ix = xStep; ix < particlesData.size().x(); ix += 4) {
                for (int iy = yStep; iy < particlesData.size().y(); iy += 4) {
                    for (int iz = 0; iz < particlesData.size().z(); ++iz) {
                        for (auto& particle : particlesData(ix,iy,iz)) {
                            const auto coord = particle.coord;
                            mesh.update_Lmat(coord, charge, _mass, _mpw);
                        }
                    }
                }
            }
#pragma omp barrier
        }
    }
}
}

