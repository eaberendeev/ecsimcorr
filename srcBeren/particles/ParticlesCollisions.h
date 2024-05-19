#ifndef PARTICLESCOLLISIONS_H_
#define PARTICLESCOLLISIONS_H_

void collision_fields( ParticlesArray& particles, ParticlesArray& particles_e, ParticlesArray& particles_i, const Mesh& mesh, const World &world, \
		    int timestep) ;
double getW_TI(double pot_I, double pot_k,double Eabs);
double getW_BM(double pot_I,double Eabs);
double getW_BS(double pot_I,double Eabs);
#endif