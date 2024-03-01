#ifndef DAMPING_H_
#define DAMPING_H_
#include "Vec.h"
#include "World.h"
void Damping_Func(double& source, double i, double maxi, double& energyDamp);
double damping_fields_circleXY(Field3d& fieldE, Field3d& fieldB,
                               const Region& domain);
double damping_fields(Field3d& fieldE, Field3d& fieldB, const Region& domain,
                      bool dampX, bool dampY, bool dampZ);
#endif 
