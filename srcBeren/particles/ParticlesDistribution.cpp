#include "Particles.h"
#include "World.h"
#include "Shape.h"

inline double pow2(double a){
  return a*a;
}
inline double pow3(double a){
  return a*a*a;
}
inline double pow4(double a){
  return a*a*a*a;
}

inline double smooth_sin(double dist){
  if (fabs(dist) > 1.) return 0;
  else
  return 0.5*(1.+cos(dist*PI));
}

double set_density_neutron(double3 r){
	return exp( -pow4( ( r.x() ) / (1.2 * Dx * PlasmaCellsX_glob / 2.) ));
//	return exp( -pow4( (x(0)-Dz*PlasmaCellsZ_glob/2)/(1.2*Dz*PlasmaCellsZ_glob/2) ));
}


void ParticlesArray::set_space_distribution(){
	
	int startX = 55 ; //0.5*_world.regionGlob.numCells.x() + 0.5*_world.regionGlob.dampCells[0].x();
	int endX = 65 ; // _world.regionGlob.numCells.x()-65 ; //0.5*_world.regionGlob.numCells.x() - _world.regionGlob.dampCells[1].x();
	int startY = 0.5*(_world.regionGlob.numCells.y() - int(widthY/Dy));
	int endY = 0.5*(_world.regionGlob.numCells.y() + int(widthY/Dy));
	int startZ = 0.5*(_world.regionGlob.numCells.z() - int(widthZ/Dz));
	int endZ = 0.5*(_world.regionGlob.numCells.z() + int(widthZ/Dz));

	if( initDist == "StrictUniform"){
		if (NumPartPerLine*NumPartPerLine*NumPartPerLine != NumPartPerCell){
			std::cout << "Wrong value NumPartPerLine! Please use NumPartPerCell = NumPartPerLine**3\n";
			exit(-1);
		}
		set_strict_uniform(int3(startX,startY,startZ),int3(endX,endY,endZ) );
	}
	if( initDist == "Uniform"){
		    SetRandSeed(200);
		set_uniform(int3(startX,startY,startZ),int3(endX,endY,endZ) );
	}

	if( initDist == "StrictUniformCircle"){
		if (NumPartPerLine*NumPartPerLine*NumPartPerLine != NumPartPerCell){
			std::cout << "Wrong value NumPartPerLine! Please use NumPartPerCell = NumPartPerLine**3\n";
			exit(-1);
		}		
		set_strict_uniform(int3(startX,startY,startZ),int3(endX,endY,endZ) );
		set_uniform_circle( int3(startX,startY,startZ),int3(endX,endY,endZ) );
	}
	if( initDist == "UniformCircle"){
		set_uniform(int3(startX,startY,startZ),int3(endX,endY,endZ) );
		set_uniform_circle( int3(startX,startY,startZ),int3(endX,endY,endZ) );
	}

	// if( initDist == "UniformCosX_dn_k"){
	// 	if (NumPartPerLine*NumPartPerLine*NumPartPerLine != NumPartPerCell){
	// 		std::cout << "Wrong value NumPartPerLine! Please use NumPartPerCell = NumPartPerLine**3\n";
	// 		exit(-1);
	// 	}
	// 	set_strict_cosX( int3(startX,startY,startZ),int3(endX,endY,endZ), distParams[0], distParams[1]);
	// }
				std::cout<< "Distribution " << name << " : "  <<initDist<<"\n";
	return;
		#ifdef PARTICLE_MASS
			if (option.smoothMass == 1){
				set_smooth_mass();
			}
			else{
					for (auto k = 0; k < size(); ++k){
						particlesData(k).mass = _mass;
					}
			}
		#endif
		
	//	if (charge == 0) 
	//		set_distribution_density(world.region,set_density_neutron);
}
void  ParticlesArray::set_pulse_distribution(){
	
	auto sigma = temperature;
	double3 pulse;
    for(auto k = 0; k < size(); ++k){
        for(auto& particle : particlesData(k)){     
		pulse.x() = Gauss(sigma / sqrt(_mass));
		pulse.y() = Gauss(sigma / sqrt(_mass));
		pulse.z() = Gauss(sigma / sqrt(_mass));
			particle.velocity = pulse; // / sqrt(_mass * _mass + dot(pulse,pulse) );;		
		} 
	}
}




void ParticlesArray::set_strict_uniform(int3 start, int3 end){
	//#if NumPartPerLine != 0
    const double3 delta = double3(Dx / NumPartPerLine,Dy / NumPartPerLine,Dz / NumPartPerLine);

    Particle particle;
    double x,y,z;
    
    x = -0.5 * delta.x();
	//std::cout << start << " "  << end << " "<< NumPartPerLine << "\n"; 
    for( auto i = 0; i < NumPartPerLine * (end.x() - start.x() ); ++i){
		x += delta.x();
		y = -0.5 * delta.y() + Dy * start.y();
		if( ! _world.region.in_region(x + Dx * start.x() ) ) continue;

		    for ( auto j = 0; j < NumPartPerLine * (end.y()-start.y() ); ++j){
				y += delta.y();
				z = -0.5 * delta.z() + Dz * start.z();
				for ( auto k = 0; k < NumPartPerLine * (end.z()-start.z() ); ++k){
					z += delta.z();

					particle.coord.x() = x + Dx * start.x() - _world.region.origin;
					particle.coord.y() = y;
					particle.coord.z() = z;
					#ifdef PARTICLE_MPW
						particle.mpw = density * Dx * Dy * Dz / (NumPartPerCell);
					#endif
					//std::cout<< particle.coord.x() << " " << particle.coord.y() <<  " " << particle.coord.z() <<  " " << particlesData.size() <<"\n";

					add_particle(particle);
				}
			}

    }
    update_count_in_cell();
   // #endif
}

void ParticlesArray::set_uniform(int3 start, int3 end){
    Particle particle;
    double x,y,z;
    
    for( auto i = 0; i < NumPartPerCell * (end.y()-start.y() ) * (end.x() - start.x() ) * (end.z()-start.z() ); ++i){
		x =  Dx * ( start.x() +  Uniform01()*(end.x()-start.x() ) );
		y =  Dy * ( start.y() +  Uniform01()*(end.y()-start.y() ) );
		z =  Dz * ( start.z() +  Uniform01()*(end.z()-start.z() ) );
		if( ! _world.region.in_region(x) ) continue;

					particle.coord.x() = x - _world.region.origin;
					particle.coord.y() = y;
					particle.coord.z() = z;
					#ifdef PARTICLE_MPW
						particle.mpw = density * Dx * Dy * Dz / (NumPartPerCell);
					#endif

					add_particle(particle);

    }
	    update_count_in_cell();

}


void ParticlesArray::set_uniform_circle(int3 start, int3 end){
    double y,z;
    double center_y = 0.5*Dy*(end.y()+start.y() );
    double center_z = 0.5*Dz*(end.z()+start.z() );
    double radius_y = 0.5 * Dy*(end.y()-start.y() ); 
    double radius_z = 0.5 * Dz*(end.z()-start.z() );    
  	
    for ( int ix = 0; ix < particlesData.size().x(); ++ix){
		for ( int iy = 0; iy < particlesData.size().y(); ++iy){
		    for ( int iz = 0; iz < particlesData.size().z(); ++iz){
		    	int ip = 0;
            	while( ip < countInCell(ix,iy,iz) ){
                    Particle particle = particlesData(ix,iy,iz)[ip];

					y = particle.coord.y() - center_y;
					z = particle.coord.z() - center_z;
					if (y*y/(radius_y*radius_y) + z*z/(radius_z*radius_z) >1. ){
						delete_particle_runtime(ix,iy,iz,ip);
					}
					else {
						ip++;
					}
				}
			}
		}
	}
    update_count_in_cell();
}


// void ParticlesArray::set_strict_cosX(int3 start, int3 end, double delta_n0, double k0){
// 	#if NumPartPerLine != 0

//     const double3 delta = double3(Dx / NumPartPerLine,Dy / NumPartPerLine,Dz / NumPartPerLine);
// 	const int Nd = 10000000;
//     const double hd = Dx * (end.x() - start.x()) / Nd;
//     double center_y = 0.5*Dy*(end.y()+start.y() );
//     double center_z = 0.5*Dz*(end.z()+start.z() );
//     double radius_y = 0.5 * Dy*(end.y()-start.y() ); 
//     double radius_z = 0.5 * Dz*(end.z()-start.z() );    
//     Particle particle;
//     double x,y,z,yy,zz;

//     double SS = 0.;
//     auto NpX = NumPartPerLine  * (end.x() - start.x());   
//     for (auto i = 0; i < Nd; ++i){
// 		SS += ( (1.0 + delta_n0 * cos( i * k0 * hd)) * hd);
//     } 
//     auto mySS = SS / NpX;
      
//     x = -0.5 * mySS;    

//     //x = -0.5 * delta.x();

//     for( auto i = 0; i < NumPartPerLine * (end.x() - start.x() ); ++i){
// 		x += (mySS / (1.0 + delta_n0 * cos( k0 * x) ) );
// 		//x += delta.x();
// 		y = -0.5 * delta.y() + Dy * start.y();
// 		if( ! _world.region.in_region(x + Dx * start.x() ) ) continue;

// 		    for ( auto j = 0; j < NumPartPerLine * (end.y()-start.y() ); ++j){
// 				y += delta.y();
// 				z = -0.5 * delta.z() + Dz * start.z();
// 				for ( auto k = 0; k < NumPartPerLine * (end.z()-start.z() ); ++k){
// 					z += delta.z();

// 					particle.coord.x() = x + Dx * start.x() - _world.region.origin;
// 					particle.coord.y() = y;
// 					particle.coord.z() = z;
// 					#ifdef PARTICLE_MPW
// 						particle.mpw = density * Dx * Dy * Dz / (NumPartPerCell);
// 					#endif
// 				//	std::cout<< particle.coord.x() << " " << particle.coord.y() <<  " " << particle.coord.z() <<  " " << particlesData.size() <<"\n";
// 					yy = particle.coord.y() - center_y;
// 					zz = particle.coord.z() - center_z;
// 					if (yy*yy/(radius_y*radius_y) + zz*zz/(radius_z*radius_z) <1. ){
// 						add_particle_scatter(particle);
// 					}
// 				}
// 			}

//     }
// #endif
// }


// void ParticlesArray::set_smooth_mass(){
// 	#ifdef PARTICLE_MASS
// 	    double x;
// 	    double3 r;

// 		for (auto k = 0; k < size(); ++k){
// 			r =  _world.region.get_coord_glob(particlesData(k).coord);  
// 			if (r.x() < 0.5 * Dx* NumCellsX_glob) 
// 				x = r.x() - Dx * DampCellsX_glob[0];
// 			else
// 				x = Dx * NumCellsX_glob - Dx * DampCellsX_glob[1] - r.x();
// 			particlesData(k).mass = (option.smoothMassMax * smooth_sin( x  / option.smoothMassSize ) + 1. ) * _mass;
// 		}
// 	#else
// 		std::cout << "Class Particle has not a mass\n";
// 	#endif
// }

// void ParticlesArray::set_distribution_density(std::function<double(double3 )> set_density){
// 	#ifdef PARTICLE_MPW
//     double x, y, z, f;
//     double3 r;

// 	for (auto k = 0; k < size(); ++k){
// 		x = particlesData(k).coord.x() - Dx*DampCellsX_glob[0];
// 		y = particlesData(k).coord.y();
// 		z = particlesData(k).coord.z();
// 		r =  _world.region.get_coord_glob(double3(x,y,z));  
// 		f = set_density(r);
// 		particlesData(k).mpw *= f;
// 	}
// 	#endif
		
// }
