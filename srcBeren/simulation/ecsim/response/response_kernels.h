#pragma once

#include "response.h"
#include "utils/bmatrix.h"

namespace ecsim_response {
// Per-particle kernels retain the dense Linear block layout and batched path.
void update_Lmat2(const Vector3R& coord, const Domain& domain, double charge, double mass, double mpw,
                  const Field3d& fieldB, double dt, Block& block);
template <int size>
void update_Lmat2(const Vector3R* coord, const Domain& domain, double charge, double mass, double mpw,
                  const Field3d& fieldB, double dt, Block& block);
}   // namespace ecsim_response
