#pragma once
#include "fields/field_boundaries.h"
#include "particles/particle_array.h"

namespace algorithmsECSIM {

template <ShapeType ST = SHAPE>
void predict_velocity(ParticlesArray& particles, const Field3d& fieldEp, const Field3d& fieldB, const double dt);
template <ShapeType ST = SHAPE>
void predict_current(const ParticlesArray& particles, const Field3d& fieldB, Field3d& fieldJ, const double dt);

Vector3R calc_JE_component(const Field3d& fieldE, const Field3d& fieldJ, const Grid& grid,
                           const FieldBoundaryConditions& bc);

}   // namespace algorithmsECSIM
