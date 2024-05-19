#include "Vec.h"
#include "World.h"
#include "Particles.h"

void ParticlesArray::set_test_particles()
{
    Particle particle;
    
    particle.coord.x() = Dx*10.01; //344/2.;
    particle.coord.y() = Dy*5.3; //40;
    particle.coord.z() = Dz*5.5; //40+27;
    
    if (name == "Electrons"){
        particle.velocity.x() = -0.02; // *Uniform01();
        particle.velocity.y() = -0.013; // *Uniform01();
        particle.velocity.z() = -0.013; //*Uniform01();
        //particle.pulse = particle.velocity*_mass/(1.-dot(particle.velocity,particle.velocity));
        add_particle(particle);
        // for(int i =0; i<30*30*1000; i++){
        //     particle.velocity.x() = 0.02*(1 - 2* Uniform01() );
        //     particle.velocity.y() = 0.023*(1 - 2* Uniform01() );
        //     particle.velocity.z() = 0; //0.023*(1 - 2* Uniform01() );
        //     particle.coord.x() = Dx*(4. + 20.*Uniform01()); //344/2.;
        //     particle.coord.y() = Dy*(4. + 20.*Uniform01()); //40;
        //     particle.coord.z() = Dz*(5. +  Uniform01()) ; //9. + 1.1*Uniform01(); //40+27;
        //     add_particle(particle);
        // }
    }
    if (name == "Ions"){

        particle.velocity.x() = -0.0025;
        particle.velocity.y() = 0.003;
        particle.velocity.z() = 0.;
        //particle.pulse = particle.velocity*_mass/(1.-dot(particle.velocity,particle.velocity));

    }
} 
