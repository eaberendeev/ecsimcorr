#pragma once

#include "coupling/esirkepov.h"
#include "particles/Particle.h"
#include "response.h"
#include "utils/bmatrix.h"

namespace ecsim_response {
// Predicted current of the particles of a cell, J = q w (I + R(b)) v / (1 + b^2)
// on the Linear Yee nodes, as predict_current. The buffer starts one node below
// the cell (deposition_origin).
using PredictedCurrent = CurrentBuffer<3>;

// Per-particle kernels retain the dense Linear block layout and batched path.
// With `current`, the predicted current of the particles is added to it, from
// the same magnetic field and rotation as the response.
void update_Lmat2(const Particle& particle, const Domain& domain, double charge, double mass, double mpw,
                  const Field3d& fieldB, double dt, Block& block, PredictedCurrent* current);
template <int size>
void update_Lmat2(const Particle* const* particles, const Domain& domain, double charge, double mass, double mpw,
                  const Field3d& fieldB, double dt, Block& block, PredictedCurrent* current);
}   // namespace ecsim_response
