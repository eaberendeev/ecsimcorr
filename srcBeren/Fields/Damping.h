#ifndef DAMPING_H_
#define DAMPING_H_
#include "Vec.h"
#include "World.h"
void Damping_Func(double& , int , int , double& );
double damping_fields(Field3d& fieldE, Field3d& fieldB,const Region& domain);
void SetEnergyDampLineNull();
#endif 
