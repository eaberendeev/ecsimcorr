#pragma once
#include "ParticlesArray.h"

namespace algorithmsECSIM {

void predict_velocity(ParticlesArray& particles, const Field3d& fieldEp,
                      const Field3d& fieldB,
                      const double dt, ShapeType type);
void predict_current(const ParticlesArray& particles, const Field3d& fieldB,
                     Field3d& fieldJ, const double dt,
                     ShapeType type);
}