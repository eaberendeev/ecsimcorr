#include "Particles.h"

void ParticlesArray::add_uniform_cilinder(int numParts, double r0, double z0, double3 c){
	Particle particle;
	auto sigma = temperature;
	double3 coord, pulse;
	double rx,ry,x,y,z;

	for(int k = 0; k < numParts; k++){
		
		do {
		rx = (1 - 2*Uniform01() ) * r0;
		ry = (1 - 2*Uniform01() ) * r0;

		} while ( rx*rx + ry*ry > r0*r0);

		

		x = c.x() + rx;
		y = c.y() + ry;
		z = c.z() + z0*(1 - 2*Uniform01() ); 
		particle.coord = double3(x,y,z);

		pulse.x() = Gauss(sigma / sqrt(_mass));
		pulse.y() = Gauss(sigma / sqrt(_mass));
		pulse.z() = Gauss(sigma / sqrt(_mass));
		particle.velocity = pulse;

		add_particle(particle);
	}
    update_count_in_cell();
}

void ParticlesArray::add_uniform_line(int numParts, double3 r, double3 sizeL){
	Particle particle;
	auto sigma = temperature;
	double3 coord, pulse;
	double x, y, z;

	for(int k = 0; k < numParts; k++){

		x = r.x() + sizeL.x()*Uniform01();
		y = r.y() + sizeL.y()*Uniform01(); 
		z = r.z() + sizeL.z()*Uniform01(); 
		particle.coord = double3(x,y,z);
		do {
			pulse.x() = Gauss(sigma / sqrt(_mass));
			pulse.y() = Gauss(sigma / sqrt(_mass));
			pulse.z() = Gauss(sigma / sqrt(_mass));
		} while (fabs(pulse.x()) > 3*sigma || fabs(pulse.y()) > 3*sigma  || fabs(pulse.z()) > 3*sigma);

		particle.velocity = pulse;

		add_particle(particle);
	}
    update_count_in_cell();
}
