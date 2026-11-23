#pragma once

#include <array>
#include <vector>

#include "fields/row_block.h"
#include "grid/domain.h"
#include "include/config.h"
#include "include/types.h"
#include "particles/particle_array.h"
#include "utils/containers.h"

namespace ecsim_response {
// Nodes of one component in the block of a cell (3 x 2 x 2).
constexpr int block_nodes = 12;
// Couplings of the block rows of one cell: three components per row and column.
constexpr int layout_entries = 3 * block_nodes * 3 * block_nodes;
// Block couplings at or below this magnitude are not added to the response.
constexpr double response_value_tolerance = 1e-16;

// Fixed pattern of the response on field unknowns: every column that one Linear
// particle couples to a row (offsets per axis [-1,1] for equal components,
// [-1,2] along the row component and [-2,1] along the column one). Columns map
// through unknown_index: periodic images join their owner, fixed components drop.
// Only rows flagged in active_rows get these columns (all rows when it is empty);
// the pattern of `united_with`, if given, is added to every row. Values are zero.
Operator response_pattern(const Domain& domain, const std::array<bool, 3>& periodic,
                          const std::vector<int>& unknown_index, const std::vector<int>& storage_index,
                          const std::vector<int8_t>& active_rows = {}, const Operator* united_with = nullptr);

// Flag the storage cells that hold particles of a species; `occupied` has the
// size of the species cell grid.
void mark_occupied_cells(const ParticlesArray& species, Array3D<int8_t>& occupied);

// Where the couplings of a cell block land in the rows of a response matrix. A
// regular row holds exactly the 123 columns of response_pattern() for its node,
// all distinct unknowns in storage order; a coupling of block row `r` (component
// c1) and block column `c` (component c2) then sits at
// position[((c1 * block_nodes + r) * 3 + c2) * block_nodes + c] from the row
// start, -1 for the couplings no particle has. Other rows (near walls and
// periodic joins, see response_layout()) are searched.
struct ResponseLayout {
    std::vector<char> regular;
    std::array<int, layout_entries> position;
};

ResponseLayout response_layout(const Operator& response, const Domain& domain, const std::vector<int>& unknown_index,
                               const std::vector<int>& storage_index);

// Add the response of the charged species to `response`, whose pattern comes
// from response_pattern() and whose layout is `layout`. Particles must already
// be grouped by cell, all species on the same cell grid.
// L = (dt/2) dJ_mid/dE_mid: the same Linear Yee weights and magnetic rotation
// as predict_current/predict_velocity are required for the discrete work balance.
// Select floor-based support and evaluate weights in logical cell coordinates;
// translate to ghost storage afterwards, so the halo does not affect weights.
// The species are summed per cell before the rows are written, so the rows of a
// cell are written once for all species.
// With `current`, the predicted current of predict_current (all species) is
// added to it in the same pass, from the same magnetic field.
template <ShapeType ST = SHAPE>
void accumulate_species(const Species& species, const Field3d& fieldB, const Domain& domain, double dt,
                        const std::vector<int>& unknown_index, const ResponseLayout& layout, Operator& response,
                        Field3d* current = nullptr);
}   // namespace ecsim_response
