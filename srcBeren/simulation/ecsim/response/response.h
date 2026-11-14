#pragma once

#include <array>
#include <vector>

#include "fields/row_block.h"
#include "grid/domain.h"
#include "include/config.h"
#include "include/types.h"
#include "utils/containers.h"

class ParticlesArray;

namespace ecsim_response {
// Fixed pattern of the response on field unknowns: every column that one Linear
// particle couples to a row (offsets per axis [-1,1] for equal components,
// [-1,2] along the row component and [-2,1] along the column one). Columns map
// through unknown_index: periodic images join their owner, fixed components drop.
Operator response_pattern(const Domain& domain, const std::array<bool, 3>& periodic,
                          const std::vector<int>& unknown_index, const std::vector<int>& storage_index);

// Add the response of a species to `response`, whose pattern comes from
// response_pattern(). Particles must already be grouped by cell; neutral
// species do not contribute.
// L = (dt/2) dJ_mid/dE_mid: the same Linear Yee weights and magnetic rotation
// as predict_current/predict_velocity are required for the discrete work balance.
// Select floor-based support and evaluate weights in logical cell coordinates;
// translate to ghost storage afterwards, so the halo does not affect weights.
template <ShapeType ST = SHAPE>
void accumulate_species(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, double dt,
                        const std::vector<int>& unknown_index, Operator& response);
}   // namespace ecsim_response
