#include "Vec.h"
#include "World.h"
#include "Particles.h"
#include "Shape.h"
#include "bounds.h"

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
    //Je(world.region.numNodes),kineticEPredict(world.region.numNodes),
    currentOnGrid(world.region.numNodes,3),
    Pxx(world.region.numNodes(0),world.region.numNodes(1)),
    Pyy(world.region.numNodes(0),world.region.numNodes(1)),
    _world(world) 
{
    for (const auto& line: vecStringParams){
        set_params_from_string(line);
    }
    injectionEnergy = 0.;


    if(RECOVERY > 0){
        read_particles_from_recovery();
        std::cout << "Upload " + name + " success!\n";
    
    } else {
        if( k_particles_reservation > 0. ){
            for(auto k = 0; k < size(); ++k){
                    particlesData(k).reserve(int(k_particles_reservation*NumPartPerCell));
            }
        }
        set_distribution();
    }
    //set_test_particles();
    // add_particle_scatter(particle);
    update_count_in_cell();
    density_on_grid_update();
    //phase_on_grid_update(); 
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
		set_pulse_distribution();	
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
    //double gama;
#pragma omp parallel for reduction(+:energy)
    for(auto k = 0; k < size(); ++k){
        for(const auto& particle : particlesData(k)){
            const auto velocity = particle.velocity;
            //pulse = velocity*_mass/(1.-dot(velocity,velocity));
            //gama = sqrt(mass(k)*mass(k) + dot(pulse,pulse) );
            energy += 0.5 * _mpw * _mass *dot(velocity,velocity); //(gama - mass(k));
        }
    }

    return energy;
}

void ParticlesArray::save_init_coord() {
//#pragma omp parallel for
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){
            particle.initCoord = particle.coord;
            particle.initVelocity = particle.velocity;
        }
    }
}

void ParticlesArray::delete_bounds(){

    Particle particle;
    for ( int ix = 0; ix < particlesData.size().x(); ++ix){
		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
		    	int ip = 0;
            	while( ip < countInCell(ix,iy,iz) ){  
			int cent = particlesData.size().x()/2;
                    if( (ix - cent)*(ix-cent)+(iy - cent)*(iy-cent) <= (cent-4)*(cent-4)){
					    ip++;
				    } else {
					    delete_particle_runtime(ix,iy,iz,ip);
                            lostEnergy += 0.5 * _mpw * _mass * dot(particle.velocity,particle.velocity); 
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


// void ParticlesArray::move(double dt){

//     Particle particle;
//     int ix2,iy2,iz2;
//     for ( int ix = 0; ix < particlesData.size().x(); ++ix){
// 		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
// 		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
// 		    	int ip = 0;
//             	while( ip < countInCell(ix,iy,iz) ){  
//                     particle = particlesData(ix,iy,iz)[ip];  
//                     particle.move(dt);
//                     ix2 = int(particle.coord.x() / Dx + CELLS_SHIFT);
//                     iy2 = int(particle.coord.y() / Dy + CELLS_SHIFT);
//                     iz2 = int(particle.coord.z() / Dz + CELLS_SHIFT);
//                     if(ix == ix2 && iy == iy2 && iz == iz2){
//                         particlesData(ix,iy,iz)[ip] = particle;  
// 					    ip++;
// 				    } else {
// 					    delete_particle_runtime(ix,iy,iz,ip);
//                         particle_boundaries(particle.coord);
//                         ix2 = int(particle.coord.x() / Dx + CELLS_SHIFT);
//                         iy2 = int(particle.coord.y() / Dy + CELLS_SHIFT);
//                         iz2 = int(particle.coord.z() / Dz + CELLS_SHIFT);
// 						particlesData(ix2,iy2,iz2).push_back(particle);
//                     }
//                 }
//             }
//         }
//     }
//     update_count_in_cell();
// }

void ParticlesArray::prepare(int timestep){
    int seed = timestep;
    //Je.clear();
    currentOnGrid.clear();
    Pxx.clear();
    Pyy.clear();
    //kineticEPredict.clear();
    lostEnergy = 0.;
    delete_bounds();

    SetRandSeed(seed%1000);
    seed+=1;

/*    double3 r0, sizeL;
    double width = Dx*125;
    r0.x() = Dx * NumCellsX_glob / 2 - 0.5*width;
    r0.y() = 0;
    r0.z() = 0;
    sizeL.x() = width;
    sizeL.y() = Dy * NumCellsY_glob;
    sizeL.z() = Dz * NumCellsZ_glob;

    add_uniform_line(125*25/2, r0, sizeL);
*/

    double3 center;
    center.x() = Dx * NumCellsX_glob / 2;
    center.y() = Dx * NumCellsY_glob / 2;
    center.z() = 0;
    add_uniform_cilinder(7065, 30*Dx, Dz * NumCellsZ_glob, center);



        int counter = 0;
        for(auto cell = 0; cell < size(); ++cell){
            for(const auto& particle : particlesData(cell)){
                counter++;   
            }
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
//#pragma omp parallel
//{
    constexpr auto SMAX = 2*SHAPE_SIZE;

    const double conx = Dx / (6*dt) * _mpw;
    const double cony = Dy / (6*dt) * _mpw;
    const double conz = Dz / (6*dt) * _mpw;

//    for ( int pStep = 0; pStep < 4; pStep++){
//#pragma omp for
        // for ( int ix = 0; ix < particlesData.size().x(); ix++){
        // for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		//     for ( int iz = 0; iz < particlesData.size().z(); ++iz){
        //     	for (auto& particle : particlesData(ix,iy,iz)){
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
        //}
        //}
//#pragma omp barrier
//    }

//} // end parallel

}

void ParticlesArray::get_P(){
    Pxx.clear();
    Pyy.clear();

// #pragma omp parallel
// {
    constexpr auto SMAX = SHAPE_SIZE;

//    for ( int pStep = 0; pStep < 4; pStep++){
//#pragma omp for
        // for ( int ix = 0; ix < particlesData.size().x(); ix++){
		// for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		//     for ( int iz = 0; iz < particlesData.size().z(); ++iz){
        //     	for (auto& particle : particlesData(ix,iy,iz)){
#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double wx[2];
        double wy[2];
        for(auto& particle : particlesData(pk)){        
                    const auto coord = particle.coord;
                    const auto velocity = particle.velocity;

                    auto x = coord.x() / Dx + CELLS_SHIFT;
                    auto y = coord.y() / Dy + CELLS_SHIFT;

                    const auto intx = int(x);
                    const auto inty = int(y);

                    wx[1] = (x - intx);
                    wx[0] = 1 - wx[1];
                    wy[1] = (y - inty);
                    wy[0] = 1 - wy[1];
                    
                    double vx =  velocity.x();
                    double vy = velocity.y();
                    double pxx = _mass*vx*vx / (particlesData.size().z() - ADD_NODES) * _mpw;
                    double pyy = _mass*vy*vy / (particlesData.size().z() - ADD_NODES) * _mpw;


                    for(int nx = 0; nx < SMAX; ++nx){
                        const int i = intx + nx;
                        for(int ny = 0; ny < SMAX; ++ny){
                            const int j = inty  + ny;
                            const auto sx = wx[nx] * wy[ny];
                            #pragma omp atomic update
                            Pxx(i,j) += sx*pxx;          
                            #pragma omp atomic update
                            Pyy(i,j) += sx*pyy;          
                        }
                    }

                }
            }
       // }
       // }
//#pragma omp barrier
//    }
//} // end parallel

}

void ParticlesArray::get_Pr(){
    Pxx.clear();
    Pyy.clear();

// #pragma omp parallel
// {
    constexpr auto SMAX = SHAPE_SIZE;

#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        double wx[2];
        double wy[2];
        for (auto& particle : particlesData(pk)){

            const auto coord = particle.coord;
            const auto velocity = particle.velocity;

            auto x = coord.x() / Dx + CELLS_SHIFT;
            auto y = coord.y() / Dy + CELLS_SHIFT;

            const auto intx = int(x);
            const auto inty = int(y);

            wx[1] = (x - intx);
            wx[0] = 1 - wx[1];
            wy[1] = (y - inty);
            wy[0] = 1 - wy[1];
                    
            double x0 = 0.5 * Dx * NumCellsX_glob;
            double R = sqrt((coord.x()-x0)*(coord.x()-x0) + (coord.y()-x0)*(coord.y()-x0) );

            double vr =   ((coord.x()-x0) / R) * velocity.x() + ((coord.y()-x0) / R) * velocity.y();
            double vp = - ((coord.y()-x0) / R) * velocity.x() + ((coord.x()-x0) / R) * velocity.y();
            double prr = _mass*vr*vr / (particlesData.size().z() - ADD_NODES) * _mpw;
            double ppp = _mass*vp*vp / (particlesData.size().z() - ADD_NODES) * _mpw;


            for(int nx = 0; nx < SMAX; ++nx){
                const int i = intx + nx;
               for(int ny = 0; ny < SMAX; ++ny){
                    const int j = inty  + ny;
                    const auto sx = wx[nx] * wy[ny];
                    #pragma omp atomic update
                    Pxx(i,j) += sx*prr;          
                    #pragma omp atomic update
                    Pyy(i,j) += sx*ppp;          
                }
            }
        }
    }

}

void ParticlesArray::move_and_calc_current(const double dt, Field3d& fieldJ){
// #pragma omp parallel
// {
    constexpr auto SMAX = 2*SHAPE_SIZE;

    const double conx = Dx / (6*dt) * _mpw;
    const double cony = Dy / (6*dt) * _mpw;
    const double conz = Dz / (6*dt) * _mpw;

//    for ( int pStep = 0; pStep < 4; pStep++){
//#pragma omp for
        // for ( int ix = 0; ix < particlesData.size().x(); ix++){
        // for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		//     for ( int iz = 0; iz < particlesData.size().z(); ++iz){
        //     	for (auto& particle : particlesData(ix,iy,iz)){
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
                    //if(fabs((end-start).z() ) >= 0.5*Dz) std::cout << particle << " " << particle.initCoord << "\n";

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
                    // ix2 = int(end.x() / Dx + CELLS_SHIFT);
                    // iy2 = int(end.y() / Dy + CELLS_SHIFT);
                    // iz2 = int(end.z() / Dz + CELLS_SHIFT);
                    // if(ix == ix2 && iy == iy2 && iz == iz2){
                    //     particlesData(ix,iy,iz)[ip] = particle;  
					//     ip++;
				    // } else {
					//     delete_particle_runtime(ix,iy,iz,ip);
                    //     particle_boundaries(particle.coord);
                    //     ix2 = int(particle.coord.x() / Dx + CELLS_SHIFT);
                    //     iy2 = int(particle.coord.y() / Dy + CELLS_SHIFT);
                    //     iz2 = int(particle.coord.z() / Dz + CELLS_SHIFT);
					// 	particlesData(ix2,iy2,iz2).push_back(particle);
                    // }
                }
            }
        //}
        //}
//#pragma omp barrier
//    }
    //update_count_in_cell();
//} // end parallel
}

void ParticlesArray::calc_Esirkepov_current(const double dt, Field3d& fieldJ) const{
    constexpr auto SMAX = 2*SHAPE_SIZE;
        const double conx = Dx / (6*dt) * _mpw;
        const double cony = Dy / (6*dt) * _mpw;
        const double conz = Dz / (6*dt) * _mpw;
    // for ( int ix = 0; ix < particlesData.size().x(); ++ix){
	// 	for ( int iy = 0; iy < particlesData.size().y(); ++iy){
	// 	    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
    //         	for (const auto& particle : particlesData(ix,iy,iz) ){
//#pragma omp parallel for
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
                    //std::cout << start << " " << end << "\n";
                    //if(fabs((end-start).z() ) >= 0.5*Dz) std::cout << particle << " " << particle.initCoord << "\n";

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
      //  }
    //}
}
// void ParticlesArray::correctv(const double dt, Field3d& fieldE, Field3d& fieldEn, Field3d& fieldEp, Field3d& fieldJp){
//     constexpr auto SMAX = 2*SHAPE_SIZE;
//     int xk, yk, zk, n, m, k,indx, indy,indz;
//     double xx, yy, zz, arg;
//     double xn, yn,zn;
//     alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
//     alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
//     alignas(64) double jx[SMAX][SMAX][SMAX];
//     alignas(64) double jy[SMAX][SMAX][SMAX];
//     alignas(64) double jz[SMAX][SMAX][SMAX];

//     const double conx = Dx / (6*dt) * density;
//     const double cony = Dy / (6*dt) * density;
//     const double conz = Dz / (6*dt) * density;
//     double3 start, end;
//     std::array<double,20> ldistr;
//     for (auto& val: ldistr){
//         val =0.0;
//     }
//     int i,j;
//     int i05,j05,k05;
//     int nx,ny,nz;
//     int ix05,iy05,iz05;

//     double q = charge;
//     int ix2,iy2,iz2;
//     for ( int ix = 0; ix < particlesData.size().x(); ++ix){
// 		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
// 		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
//             	for (auto& particle : particlesData(ix,iy,iz) ){
//                     start = particle.initCoord;                    
//                     end = particle.coord;

//                     xx = start.x() / Dx;
//                     yy = start.y() / Dy;
//                     zz = start.z() / Dz;

//                     xn = end.x() / Dx;
//                     yn = end.y() / Dy;
//                     zn = end.z() / Dz;

//                     xk = int(xx);
//                     yk = int(yy);
//                     zk = int(zz);

//                     for(n = 0; n < SMAX; ++n){
//                         for(m = 0; m < SMAX; ++m){
//                             for(k = 0; k < SMAX; ++k){
//                                 jx[n][m][k] = 0.;
//                                 jy[n][m][k] = 0.;
//                                 jz[n][m][k] = 0.;
//                             }
//                         }
//                     }

//                     for(n = 0; n < SMAX; ++n){
//                         arg = -xx + double(xk - CELLS_SHIFT + n);
//                         sx[n] = Shape(arg);
//                         arg = -yy + double(yk - CELLS_SHIFT + n);
//                         sy[n] = Shape(arg);
//                         arg = -zz + double(zk - CELLS_SHIFT + n);
//                         sz[n] = Shape(arg);            
//                         arg = -xn + double(xk - CELLS_SHIFT + n);
//                         sx_n[n] = Shape(arg);
//                         arg = -yn + double(yk - CELLS_SHIFT + n);
//                         sy_n[n] = Shape(arg);
//                         arg = -zn + double(zk - CELLS_SHIFT + n);
//                         sz_n[n] = Shape(arg);
//                     }

//                     for(n = 0; n < SMAX; ++n){
//                         indx = xk  + n;
//                         for(m = 0; m < SMAX; ++m){
//                             indy = yk + m;
//                             for(k = 0; k < SMAX; ++k){
//                                 indz = zk + k;
                        
//                                 if(n == 0) jx[n][m][k] = -q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

//                                 if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
                                
//                                 if(m == 0) jy[n][m][k] = -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
//                                 if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
                                
//                                 if(k == 0) jz[n][m][k] = -q * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
//                                 if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-q * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
//                             }
//                         }
//                     }
//                     double je_cell(0);
//                     for(n = 0; n < SMAX; ++n){
//                         indx = xk  + n;
//                         for(m = 0; m < SMAX; ++m){
//                             indy = yk + m;
//                             for(k = 0; k < SMAX; ++k){
//                                 indz = zk + k;
//                                 je_cell += 0.5*(fieldEp(indx,indy,indz,0) + fieldEn(indx,indy,indz,0)) * jx[n][m][k] +
//                                     0.5*(fieldEp(indx,indy,indz,1) + fieldEn(indx,indy,indz,1)) * jy[n][m][k] +
//                                     0.5*(fieldEp(indx,indy,indz,2) + fieldEn(indx,indy,indz,2)) * jz[n][m][k];
//                             }
//                         }
//                     }
//                     const auto initVelocity = particle.initVelocity;
//                     const auto velocity = particle.velocity;

//                     double3 coord = start + 0.5*Dt*initVelocity;
                    
//                     double3 Ep = get_fieldE_in_pos(fieldEp,coord); 
//                     double3 En = get_fieldE_in_pos(fieldEn,coord);

//                     double3 v12 = 0.5*(velocity + initVelocity); 

//                     double jp_cell(0);
//                     jp_cell = 0.5*charge*dot(v12,(Ep+En));

//                    // std::cout << je_cell << " " << jp_cell << "\n";

//                     double energyPredict = 0.5 * _mass *dot(velocity,velocity);

//                     Je(ix,iy,iz) += Dt*(je_cell - jp_cell);
                    
//                     kineticEPredict(ix,iy,iz) += energyPredict;

//                 }
//             }
//         }
//     }
//     for ( int ix = 0; ix < particlesData.size().x(); ++ix){
// 		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
// 		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
//             	for (auto& particle : particlesData(ix,iy,iz) ){

//                     const auto initVelocity = particle.initVelocity;
//                     const auto velocity = particle.velocity;
//                     const double energyPredict = kineticEPredict(ix,iy,iz);
//                     const double lambda = sqrt(1 + Je(ix,iy,iz) / energyPredict);
//                     // 0. - 2. h = 0.1
//                     int indL = int(lambda / 0.1);
//                     ldistr[indL]+=1.;
//                     if(fabs(1-lambda) > 0.01)
//                     std::cout << "lambda "<< lambda << " v_start " << initVelocity 
//                     << " v_fin " <<  lambda*velocity << " v_predict " <<  particle.velocity 
//                     << " 0.5*mv^2_predict " << energyPredict << " JE " << Je(ix,iy,iz) << "\n";
//                     particle.velocity = lambda*velocity;
//                 }
//             }
//         }
//     }

// }


void ParticlesArray::correctv(Mesh& mesh){

    std::array<double,20> ldistr;
    for (auto& val: ldistr){
        val =0.0;
    }

    double jp_cell = 0;
#pragma omp parallel for reduction(+:jp_cell)
    // for ( int ix = 0; ix < particlesData.size().x(); ++ix){
	// 	for ( int iy = 0; iy < particlesData.size().y(); ++iy){
	// 	    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
    //         	for (auto& particle : particlesData(ix,iy,iz) ){
//#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    const auto initVelocity = particle.initVelocity;
                    const auto velocity = particle.velocity;
 
                    double3 start = particle.initCoord;
                    double3 coord = start + 0.5*Dt*initVelocity;
                    
                    double3 Ep = get_fieldE_in_pos(mesh.fieldEp,coord); 
                    double3 E = get_fieldE_in_pos(mesh.fieldE,coord);

                    double3 v12 = 0.5*(velocity + initVelocity); 

                    jp_cell += 0.5*_mpw*charge*dot(v12,(Ep+E));
                    //std::cout << "jp " << jp_cell  << " " << Ep << "\n";

                }
            }
      //  }
    //}
    
    
    mesh.fieldJp_full.data() = mesh.fieldJp.data() + mesh.Lmat2*(mesh.fieldE.data() + mesh.fieldEp.data())/Dt;
    //const double energyEn = mesh.calc_JE(mesh.fieldEn,mesh.fieldJe);
    //const double energyEp = mesh.calc_JE(mesh.fieldEp,mesh.fieldJp_full);

    // mesh.fieldJp_full.data() = mesh.fieldJp.data() + mesh.Lmat2*(mesh.fieldE.data() + mesh.fieldEp.data())/Dt;
    const double energyJeEn = mesh.calc_JE(mesh.fieldEn,currentOnGrid);
    const double energyJeE = mesh.calc_JE(mesh.fieldE,currentOnGrid);
    const double energyJpEp = mesh.calc_JE(mesh.fieldEp,mesh.fieldJp_full);
    const double energyJpE = mesh.calc_JE(mesh.fieldE,mesh.fieldJp_full);
    //std::cout << "Je vs Lmat " << jp_cell << " " << 0.5*(energyJpE + energyJpEp) << "\n";
    const double energyK = get_kinetic_energy();
    const double lambda = sqrt(1 + Dt*( 0.5*(energyJeEn + energyJeE) - jp_cell) / energyK);

// #pragma omp parallel for
//     for ( int ix = 0; ix < particlesData.size().x(); ++ix){
// 		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
// 		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
//             	for (auto& particle : particlesData(ix,iy,iz) ){
#pragma omp parallel for
    for(auto pk = 0; pk < size(); ++pk){
        for(auto& particle : particlesData(pk)){     
                    //const auto initVelocity = particle.initVelocity;
                    const auto velocity = particle.velocity;

                    // if(fabs(1-lambda) > 0.005)
                    // std::cout << "lambda "<< lambda << " v_start " << initVelocity 
                    // << " v_fin " <<  lambda*velocity << " v_predict " <<  particle.velocity 
                    // //<< " 0.5*mv^2_predict " << energyPredict << " JE " << Je(ix,iy,iz) 
                    // << "\n";
                    particle.velocity = lambda*velocity;
                }
            }
      //  }
    //}
                        // 0. - 2. h = 0.1
                    //int indL = int(lambda / 0.1);
                    //ldistr[indL]+=1.;
    const double energyK2 = get_kinetic_energy();
    std::cout << "lambda "<< lambda  << " " << lambda*lambda << " "<< energyK2-energyK << " " << 0.5*Dt*(energyJeEn + energyJeE - energyJpEp - energyJpE) << "\n";
   // std::cout << lambda*lambda*energyK << " " <<  energyK - (energyEn - energyEp)<< "\n";

}

void ParticlesArray::predict_current(const Field3d& fieldB, Field3d& fieldJ){
//#pragma omp parallel
//{
    constexpr auto SMAX = SHAPE_SIZE;

    

//    for ( int pStep = 0; pStep < 4; pStep++){
//#pragma omp for
        // for ( int indx = 0; indx < particlesData.size().x(); indx++){
        // for ( int indy = 0; indy < particlesData.size().y(); ++indy){
		//     for ( int indz = 0; indz < particlesData.size().z(); ++indz){
        //     	for (auto& particle : particlesData(indx,indy,indz)){   
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
        //}
       // }
//#pragma omp barrier
//    }
//} // end parallel
}

void ParticlesArray::correctv2(const Field3d& fieldB,const Field3d& fieldE){
    constexpr auto SMAX = SHAPE_SIZE;
    double3 h, B, coord;
    double3 velocity,current;
    int ix,iy,iz;
    int ix05,iy05,iz05;
    int i,j,k;
    int i05,j05,k05;
    int nx,ny,nz;
    double sx,sy,sz;
    double alpha,alpha2,beta;
    double x,y,z;
    double x05,y05,z05;
    double qp = charge;// / double(NumPartPerCell);

    alignas(64) double wx[SMAX], wy[SMAX], wz[SMAX];
    alignas(64) double wx05[SMAX], wy05[SMAX], wz05[SMAX];
    
    for(auto cell = 0; cell < size(); ++cell){
        for(auto& particle : particlesData(cell)){        
            coord = particle.initCoord + 0.5*Dt*particle.initVelocity;
            velocity = particle.initVelocity;
            std::cout << "cv2 " << coord << " " << velocity << "\n";
            x = coord.x() / Dx + CELLS_SHIFT;
            y = coord.y() / Dy + CELLS_SHIFT;
            z = coord.z() / Dz + CELLS_SHIFT;
            x05 = x - 0.5;
            y05 = y - 0.5;
            z05 = z - 0.5;

            ix = int(x);
            iy = int(y);
            iz = int(z);
            ix05 = int(x05);
            iy05 = int(y05);
            iz05 = int(z05);
            
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
            
            B = get_fieldB_in_pos(fieldB,coord); 

            beta = Dt * qp / _mass;
            alpha = 0.5*beta*mag(B);
            alpha2 = alpha*alpha;
            h = unit(B);

            current = qp * _mpw / (1.+alpha2) * (velocity + alpha*cross(velocity,h) + alpha2*dot(h,velocity)*h );
            std::cout << " current2 " << current << "\n";
            double je = 0;
            for(nx = 0; nx < SMAX; ++nx){
                i = ix + nx;
                i05 = ix05 + nx;
                for(ny = 0; ny < SMAX; ++ny){
                    j = iy  + ny;
                    j05 = iy05  + ny;
                    for(nz = 0; nz < SMAX; ++nz){
                        k = iz  + nz;
                        k05 = iz05  + nz;
                        sx = wx05[nx] * wy[ny] * wz[nz];
                        sy = wx[nx] * wy05[ny] * wz[nz];
                        sz = wx[nx] * wy[ny] * wz05[nz];
                        je += fieldE(i05,j,k,0) * sx*current.x() +
                        fieldE(i,j05,k,1) * sy*current.y() +
                        fieldE(i,j,k05,2) * sz*current.z();                   
                    }
                }
            }

                    const auto velocity = particle.velocity;
                    const auto initVelocity = particle.initVelocity;
                    double energyStart = 0.5*Dx*Dy*Dz * _mpw * _mass *dot(initVelocity,initVelocity);
                    energyStart += Dt*je;
                    const double energyEnd = 0.5*Dx*Dy*Dz * _mpw * _mass *dot(velocity,velocity);
                    const double lambda = sqrt(energyStart/energyEnd);
                    particle.velocity = lambda*velocity;
                    std::cout << lambda << " " << initVelocity << " " << velocity << " " <<  particle.velocity << " " << energyStart << " " << energyEnd << "\n";

        }
    }
}


// void ParticlesArray::predict_current2(const Field3d& fieldB, Field3d& fieldJ){
//     int jmax = size();
//     constexpr auto SMAX = 2*SHAPE_SIZE;

//     double3 h, B, POS;
//     double xx,yy,zz;
//     int xk,yk,zk,indx,indy,indz;
//     int m,n,k;
//     double arg;
//     double snm1,snm2,snm3,snm12,snm13,snm23;
//     double alpha,alpha2,beta;
//     double3 velocity,J;

//     alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
//     alignas(64) double sdx[SMAX], sdy[SMAX], sdz[SMAX];
//     for (int j = 0; j < jmax; j++ ) {
//         B = 0.;
//         POS = particlesData(j).coord;
//         velocity = particlesData(j).velocity;

//         xx = POS.x() / Dx;
//         yy = POS.y() / Dy;
//         zz = POS.z() / Dz;

//         xk = int(xx);
//         yk = int(yy);
//         zk = int(zz);
    
//         for(n = 0; n < SMAX; ++n){
//             arg = -xx + double(xk - CELLS_SHIFT + n);
//             sx[n] = Shape(arg);
//             sdx[n] = Shape(arg + 0.5);
//             arg = -yy + double(yk - CELLS_SHIFT + n);
//             sy[n] = Shape(arg);
//             sdy[n] = Shape(arg + 0.5);
//             arg = -zz + double(zk - CELLS_SHIFT + n);
//             sz[n] = Shape(arg);
//             sdz[n] = Shape(arg + 0.5);
//         }
        
//         for(n = 0; n < SMAX; ++n){
//                 indx = xk + n;

//             for(m = 0; m < SMAX; ++m){
//                 indy = yk  + m;
//                 for(k = 0; k < SMAX; ++k){

//                     snm12 = sdx[n] * sdy[m] * sz[k];
//                     snm13 = sdx[n] * sy[m] * sdz[k];
//                     snm23 = sx[n] * sdy[m] * sdz[k];
//                     indz = zk  + k;
//                     B.x() += (snm23 * fieldB(indx,indy,indz,0) );
//                     B.y() += (snm13 * fieldB(indx,indy,indz,1) );
//                     B.z() += (snm12 * fieldB(indx,indy,indz,2) );
//                 }
//             }
//         }
        
//         beta = Dt * charge / _mass;
//         alpha = 0.5*beta*mag(B);
//         alpha2 = alpha*alpha;
//         h = B / mag(B);

//         J = charge / (1.+alpha2) * (velocity+alpha2*dot(h,velocity)*h );
        
//         for(n = 0; n < SMAX; ++n){
//             indx = xk  + n;
//             for(m = 0; m < SMAX; ++m){
//                 indy = yk + m;
//                 for(k = 0; k < SMAX; ++k){
                
//                     indz = zk + k;

//                     snm1 = sdx[n] * sy[m] * sz[k];
//                     snm2 = sx[n] * sdy[m] * sz[k];
//                     snm3 = sx[n] * sy[m] * sdz[k];
//                     fieldJ(indx,indy,indz,0) += snm1*J.x();
//                     fieldJ(indx,indy,indz,1) += snm2*J.y();
//                     fieldJ(indx,indy,indz,2) += snm3*J.z();
//                 }
//             }
//         }

//     }
// }


// void ParticlesArray::stencil_Lmat( Mesh& mesh)
// {
//         const int SMAX = SHAPE_SIZE;
//     std::vector<Trip> trips;
//     const auto sizeM = mesh.fieldE.size();
//     int totalSize = sizeM.x()*sizeM.y()*sizeM.z()*9*SHAPE_SIZE*SHAPE_SIZE*SHAPE_SIZE;
//     trips.reserve(totalSize);
//     int sizeParticles = size();
//     double3 B,h;
//     double3 coord;
//     double alpha;
//     double Ap ;
//     double wx,wy,wz,wx1,wy1,wz1;
//     double value; //,column;
//     int cellLocX,cellLocY,cellLocZ,cellLocX05,cellLocY05,cellLocZ05;
//     double coordLocX,coordLocY,coordLocZ;
//     double coordLocX05,coordLocY05,coordLocZ05;
//     int i,j,k,i1,j1,k1;
//     int indx,indy,indz;
//     int vindx,vindy,vindz;
//     int indx05,indy05,indz05;
//     alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
//     alignas(64) double sx05[SMAX], sy05[SMAX], sz05[SMAX];
//     double qp = charge; // / double(NumPartPerCell);


//     for (int n = 0; n < sizeParticles; n++ ) {

//         coord = particlesData(n).coord;
//         B = 0.;

//         coordLocX = coord.x() / Dx + CELLS_SHIFT;
//         coordLocY = coord.y() / Dy + CELLS_SHIFT;
//         coordLocZ = coord.z() / Dz + CELLS_SHIFT;
//         coordLocX05 = coord.x() / Dx + CELLS_SHIFT - 0.5;
//         coordLocY05 = coord.y() / Dy + CELLS_SHIFT - 0.5;
//         coordLocZ05 = coord.z() / Dz + CELLS_SHIFT - 0.5;

//         cellLocX = int(coordLocX);
//         cellLocY = int(coordLocY);
//         cellLocZ = int(coordLocZ);
//         cellLocX05 = int(coordLocX05);
//         cellLocY05 = int(coordLocY05);
//         cellLocZ05 = int(coordLocZ05);
        
//         sx[1] = (coordLocX - cellLocX);
//         sx[0] = 1 - sx[1];
//         sy[1] = (coordLocY - cellLocY);
//         sy[0] = 1 - sy[1];
//         sz[1] = (coordLocZ - cellLocZ);
//         sz[0] = 1 - sz[1];

//         sx05[1] = (coordLocX05 - cellLocX05);
//         sx05[0] = 1 - sx05[1];
//         sy05[1] = (coordLocY05 - cellLocY05);
//         sy05[0] = 1 - sy05[1];
//         sz05[1] = (coordLocZ05 - cellLocZ05);
//         sz05[0] = 1 - sz05[1];        
        
//         for(i = 0; i < SMAX; ++i){
//             indx = cellLocX + i;
//             indx05 = cellLocX05 + i;
//             for(j = 0; j < SMAX; ++j){
//                 indy = cellLocY  + j;
//                 indy05 = cellLocY05  + j;
//                 for(k = 0; k < SMAX; ++k){
//                     indz = cellLocZ  + k;
//                     indz05 = cellLocZ05  + k;
//                     wx = sx[i] * sy05[j] * sz05[k];
//                     wy = sx05[i] * sy[j] * sz05[k];
//                     wz = sx05[i] * sy05[j] * sz[k];
//                     B.x() += wx * mesh.fieldB(indx,indy05,indz05,0);
//                     B.y() += wy * mesh.fieldB(indx05,indy,indz05,1);
//                     B.z() += wz * mesh.fieldB(indx05,indy05,indz,2);
//                 }
//             }
//         }
//         h = unit(B);
//                 // std::cout<< "B from L2 " << B << "\n";
//                 // std::cout<< "h from L2 " << h << "\n";
//         alpha = 0.5*Dt*qp*mag(B) / _mass;
//         Ap = 0.25*Dt*Dt*qp*qp / _mass / (1+alpha*alpha); 
        
//         for(i = 0; i < SMAX; ++i){
//             indx = cellLocX + i;
//             indx05 = cellLocX05 + i;
//             for(j = 0; j < SMAX; ++j){
//                 indy = cellLocY  + j;
//                 indy05 = cellLocY05  + j;
//                 for(k = 0; k < SMAX; ++k){
//                     indz = cellLocZ  + k;
//                     indz05 = cellLocZ05  + k;
//                     wx = sx05[i] * sy[j] * sz[k];
//                     wy = sx[i] * sy05[j] * sz[k];
//                     wz = sx[i] * sy[j] * sz05[k];
//                     vindx = mesh.vind(cellLocX05 + i,cellLocY + j,cellLocZ+k,0);
//                     vindy = mesh.vind(cellLocX + i,cellLocY05 + j,cellLocZ+k,1);
//                     vindz = mesh.vind(cellLocX + i,cellLocY + j,cellLocZ05+k,2);
                    
//                     for(i1 = 0; i1 < SMAX; ++i1){
//                         for(j1 = 0; j1 < SMAX; ++j1){
//                             for(k1 = 0; k1 < SMAX; ++k1){
//                                 wx1 = sx05[i1] * sy[j1] * sz[k1];
//                                 wy1 = sx[i1] * sy05[j1] * sz[k1];
//                                 wz1 = sx[i1] * sy[j1] * sz05[k1];
//                                 // xx
//                                 value = wx*wx1*Ap*(1.+alpha*alpha*h.x()*h.x() );
//                                 if(fabs(value) > 1.e-16){
//                                  //std::cout << "L2ind " << indx05 << " " << indy << " " << indz << " " << mesh.vind(indx05,indy,indz,0) <<"\n";
//                                  //std::cout << "L2ind " << cellLocX05 + i1 << " " << cellLocY + j1 << " " << cellLocZ+k1 << " " << mesh.vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0) <<"\n";                                   
//                                     mesh.Lmat2.coeffRef(vindx,mesh.vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0) ) += value;
//                                     //trips.push_back(Trip(vindx,mesh.vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0), value ));
//                                 }     
//                                 // xy
//                                 value = wx*wy1*Ap*(alpha*h.z() + alpha*alpha*h.x()*h.y() );
//                                 //std::cout << alpha << " " << wx << " " << wy1 << " " << value << "\n";
//                                 if(fabs(value) > 1.e-16){
//                                     mesh.Lmat2.coeffRef(vindx,mesh.vind(cellLocX + i1, cellLocY05 + j1, cellLocZ+k1,1) )+= value;
//                                     //std::cout << "L2ind " << vindx << " " << mesh.vind(cellLocX + i1, cellLocY05 + j1, cellLocZ+k1,1) << " " << value <<"\n";
// //                                    std::cout << "L2ind " << indx05 << " " << indy << " " << indz << " " << mesh.vind(indx05,indy,indz,0) << " " << value <<"\n";
//   //                                  std::cout << "L2ind " << cellLocX + i1 << " " << cellLocY05 + j1 << " " << cellLocZ+k1 << " " << mesh.vind(cellLocX + i1,cellLocY05 + j1,cellLocZ+k1,1) <<"\n";
//                                 }
//                                 // xz
//                                 value = wx*wz1*Ap*alpha*(-h.y() + alpha*h.x() *h.z() );
//                                  if(fabs(value) > 1.e-16){
//                                      mesh.Lmat2.coeffRef(vindx,mesh.vind(cellLocX + i1,cellLocY + j1,cellLocZ05 + k1,2) ) += value;
//                                  }
//                                 // yx
//                                 value = wy*wx1*Ap*alpha*(-h.z() + alpha*h.x() *h.y() );
//                                 if(fabs(value) > 1.e-16){
//                                     mesh.Lmat2.coeffRef(vindy,mesh.vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0)) += value;
//                                 }
//                                 // yy 
//                                 value = wy*wy1*Ap*(1.+alpha*alpha*h.y()*h.y() );
//                                 if(fabs(value) > 1.e-16){
//                                     mesh.Lmat2.coeffRef(vindy,mesh.vind(cellLocX + i1,cellLocY05 + j1,cellLocZ+k1,1)) += value;
//                                 }
//                                 // yz
//                                 value = wy*wz1*Ap*alpha*(h.x() + alpha*h.y() *h.z() );
//                                 if(fabs(value) > 1.e-16){
//                                     mesh.Lmat2.coeffRef(vindy,mesh.vind(cellLocX + i1,cellLocY + j1,cellLocZ05+k1,2)) += value;
//                                 }
//                                 // zx
//                                  value = wz*wx1*Ap*alpha*(h.y() + alpha*h.x() *h.z() );
//                                 if(fabs(value) > 1.e-16){
//                                     mesh.Lmat2.coeffRef(vindz,mesh.vind(cellLocX05 + i1,cellLocY + j1,cellLocZ+k1,0)) += value;
//                                 }
//                                 // zy
//                                 value = wz*wy1*Ap*alpha*(-h.x() + alpha*h.y() *h.z() );
//                                 if(fabs(value) > 1.e-16){
//                                     mesh.Lmat2.coeffRef(vindz,mesh.vind(cellLocX + i1,cellLocY05 + j1,cellLocZ+k1,1)) += value;
//                                 }
//                                 // zz
//                                 value = wz*wz1*Ap*(1.+alpha*alpha*h.z()*h.z() );
//                                 if(fabs(value) > 1.e-16){
//                                     mesh.Lmat2.coeffRef(vindz,mesh.vind(cellLocX + i1,cellLocY + j1,cellLocZ05+k1,2)) += value;
//                                 }
//                             } // G'z
//                         } // G'y
//                     } // G'x
//                 } // Gz
//             } // Gy
//         } // Gx
//     } /// end particles cycle

//    // mesh.Lmat2.setFromTriplets(trips.begin(), trips.end());
// }

void reserve_Lmat(Mesh& mesh, ParticlesArray &sp){
        // for ( int ix = 1; ix < particlesData.size().x()-1; ix++){
        // for ( int iy = 1; iy < particlesData.size().y()-1; ++iy){
		//     for ( int iz = 1; iz < particlesData.size().z()-1; ++iz){
        //   if(particlesData(ix,iy,iz) > 0){
#pragma omp parallel for
        for(auto pk = 0; pk < sp.size(); ++pk){
            if(sp.particlesData(pk).size() > 0){    
              if (mesh.done[pk] != 1){
                auto ix = pos_ind(pk,0, sp.particlesData.size().x(),sp.particlesData.size().y(),sp.particlesData.size().z());
                auto iy = pos_ind(pk,1, sp.particlesData.size().x(),sp.particlesData.size().y(),sp.particlesData.size().z());
                auto iz = pos_ind(pk,2, sp.particlesData.size().x(),sp.particlesData.size().y(),sp.particlesData.size().z());
                for(int ix1 = 0; ix1 < 5; ix1++){
                for(int iy1 = 0; iy1 < 5; iy1++){
                for(int iz1 = 0; iz1 < 5; iz1++){
                    double coordx = ix*Dx - 0.25*Dx+ix1*Dx;
                    double coordy = iy*Dy- 0.25*Dy+iy1*Dy;
                    double coordz = iz*Dz- 0.25*Dz+iz1*Dz;
                    mesh.reserve_Lmat(double3(coordx,coordy,coordz));
                }
                }
                }
                mesh.done[pk] = 1;
              }
          }
        }
}


void ParticlesArray::get_L(Mesh& mesh) {
#pragma omp parallel
{
    for (int xStep = 0; xStep < 4; xStep++) {
        for (int yStep = 0; yStep < 4; yStep++) {
            #pragma omp for collapse(2) schedule(dynamic)
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

