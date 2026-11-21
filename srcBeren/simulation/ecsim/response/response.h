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

// Where the couplings of a cell block land in the rows of a response matrix. A
// regular row holds exactly the 123 columns of response_pattern() for its node,
// all distinct unknowns in storage order; a coupling of block row `r` (component
// c1) and block column `c` (component c2) then sits at
// position[((c1 * 12 + r) * 3 + c2) * 12 + c] from the row start, -1 for the
// couplings no particle has. Other rows (near walls, periodic joins, I-M columns
// outside the pattern) are searched.
struct ResponseLayout {
    std::vector<char> regular;
    std::array<int, 3 * 12 * 3 * 12> position;
};

ResponseLayout response_layout(const Operator& response, const Domain& domain, const std::vector<int>& unknown_index,
                               const std::vector<int>& storage_index);

// Add the response of a species to `response`, whose pattern comes from
// response_pattern() and whose layout is `layout`. Particles must already be
// grouped by cell; neutral species do not contribute.
// L = (dt/2) dJ_mid/dE_mid: the same Linear Yee weights and magnetic rotation
// as predict_current/predict_velocity are required for the discrete work balance.
// Select floor-based support and evaluate weights in logical cell coordinates;
// translate to ghost storage afterwards, so the halo does not affect weights.
template <ShapeType ST = SHAPE>
void accumulate_species(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, double dt,
                        const std::vector<int>& unknown_index, const ResponseLayout& layout, Operator& response);
}   // namespace ecsim_response
