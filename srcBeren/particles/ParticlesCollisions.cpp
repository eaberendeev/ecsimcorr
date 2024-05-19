#include "Particles.h"
#include "Mesh.h"
#include "Diagnostic.h"
#include <cmath>
#include <algorithm> 
#include <string>
#include "ParticlesCollisions.h"
#include "Particles.h"
#include "Shape.h"

// void collision(const Mesh &mesh, const World& world ,std::vector<ParticlesArray> &Particles,int timestep){
  
//   if (!IONIZATION ) return;

//   auto num_neutrals = get_num_of_type_particles(Particles, "Neutrals");
//   auto num_e = get_num_of_type_particles(Particles, "Electrons");
//   auto num_i = get_num_of_type_particles(Particles, "Ions");
//   auto num_i2 = get_num_of_type_particles(Particles, "Ions2");
//   if (num_e >= 0 && num_i >= 0 && num_i2 >= 0)
//   	collision_fields( Particles[num_i], Particles[num_e], Particles[num_i2],mesh, world, timestep);      
//   if (num_e >= 0 && num_i >= 0 && num_neutrals >= 0)
//   	collision_fields( Particles[num_neutrals], Particles[num_e], Particles[num_i],mesh,  world, timestep);

// }

// /// Ionization particles = electron(particles_e) +  ion (particles_i) 
// void collision_fields( ParticlesArray& particlesBasic, ParticlesArray& particlesE, ParticlesArray& particlesI,
// 	                   const Mesh &mesh, const World &world, int timestep){


// //	double2 coord;
// 	double p,phi,psi,pxy;
// 	Particle particle;
// 	double3 coordGlob;
// 	double3 E_loc, E_laser;
// 	double W, W_TI, W_BM, W_BS, Prob;
// 	double Eabs;
// 	int kmax = particlesBasic.size();
// 	double Pcoll;

// 	if( kmax <= 0 ) return;
// 	int k = 0;

// 	while( k <  particlesBasic.size() ) {
// 		E_loc = get_fieldE_in_pos(mesh.fieldE,particlesBasic(k).coord);
	
// 		E_laser = 0.;
// 		// for (const auto &las: mesh.lasers){
// 		// 	coordGlob = world.region.get_coord_glob(particlesBasic(k).coord);
// 		// 	if(las.type == "Virtual"){
// 		//  		E_laser += las.get_field_coll(coordGlob, timestep);
// 		//  	}
// 		// }
		
// 		E_loc += E_laser;
// 		Eabs = mag(E_loc);

// 		W_BM = getW_BM(particlesBasic.pot_I, Eabs);
// 		W_BS = getW_BS(particlesBasic.pot_I,Eabs);

// 		W_TI = getW_TI(particlesBasic.pot_I,particlesBasic.pot_k,Eabs);

// 		W = std::min(W_TI, std::min(W_BM, W_BS));

// 		Prob = Uniform01();
// 		Pcoll = 1. - exp(-W*Dt);

// 		if( Prob <= Pcoll){
// 		    particle = particlesBasic(k);
// 		    particlesI.particlesData.push_back(particle);

// 		    p = Uniform01()*PulseFromKev(particlesBasic.pot_I, particlesE.mass() );
// 		    phi = 2.*PI*Uniform01();
// 		    psi = 2.*PI*Uniform01();
// 		    pxy = p * cos(psi);
// 		    particle.pulse.x() = pxy*sin(phi);
// 		    particle.pulse.y() = pxy*cos(phi);
// 		    particle.pulse.z() = p*sin(psi);
// 		    #ifdef PARTICLE_MASS
// 		    particle.mass = particlesE.mass();
// 		    #endif
// 		   particlesE.particlesData.push_back(particle);
// 		   particlesBasic.particlesData.del(k);

// 		}  else ++k;
// 	}

// }




double getW_TI(double pot_I, double pot_k,double Eabs){
    if (Eabs == 0.0) return 0;

  double k = sqrt(pot_I / 0.0136);
  double n = pot_k / k;
  const auto wa = 7.328e11 / sqrt(n0);
  const auto Ea = 5.356e9 / sqrt(n0);
  double F = Eabs / (k*k*k*Ea);
  double C1 = pow(2.,2.*n) / (n * std::tgamma(2.*n));
  
  
  return wa * k* k * C1 * pow(2./F,2*n-1) * exp(-2./(3*F));
}

double getW_BM(double pot_I,double Eabs){
  double k = sqrt(pot_I / 0.0136);
  const auto wa = 7.328e11 / sqrt(n0);
  const auto Ea = 5.356e9 / sqrt(n0);
  return 2.4 * wa / (k * k * k * k) * (Eabs / Ea) * (Eabs / Ea) ;
  
}

double getW_BS(double pot_I,double Eabs){
  double k = sqrt(pot_I / 0.0136);
  const auto wa = 7.328e11 / sqrt(n0);
  const auto Ea = 5.356e9 / sqrt(n0);
  
  return 0.8 * wa / k * (Eabs / Ea)  ;
  
}
