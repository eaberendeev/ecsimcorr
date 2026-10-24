#pragma once
#include "ParticlesArray.h"
#include "boundary_conditions.h"

namespace algorithmsECSIM {

template <ShapeType ST = SHAPE>
void predict_velocity(ParticlesArray& particles, const Field3d& fieldEp, const Field3d& fieldB, const double dt);
template <ShapeType ST = SHAPE>
void predict_current(const ParticlesArray& particles, const Field3d& fieldB, Field3d& fieldJ, const double dt);
void calculate_current(const ParticlesArray& particles, Field3d& fieldJ);

Vector3R calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ, const Grid& grid,
                           const BoundaryConditionHandler& bc);

}   // namespace algorithmsECSIM
