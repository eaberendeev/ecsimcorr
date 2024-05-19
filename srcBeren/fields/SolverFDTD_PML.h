#ifndef SolverFDTD_PML_H_
#define SolverFDTD_PML_H_
#include "Vec.h"
#include "World.h"
void solver_FDTD_PML(Field3d& fieldE, Field3d& fieldB,
                     Field3d& fieldEp, Field3d& fieldBp, 
                     const Field3d& fieldJ, const World& world);

#endif 	
