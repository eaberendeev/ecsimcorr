#pragma once

#include "fields/row_block.h"
#include "grid/domain.h"
#include "include/config.h"
#include "include/types.h"
#include "utils/containers.h"
#include "utils/memory.h"

class ParticlesArray;

namespace ecsim_response {
// Per-thread response contributions are assembly results, passed explicitly
// from species accumulation to CSR conversion. Reset once before all species.
using RowBlocks = SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>;
// Clear the previous assembly, adapting the buffers only if the thread count changed.
void clear_row_blocks(RowBlocks& row_blocks);

// Temporary memory used entirely within assemble_matrix(). Simulation keeps
// allocations between calls; no physical response data is carried between stages.
// All entries are overwritten on each call and are not saved in checkpoints.
struct Workspace {
    SimpleArrayBuffer<std::array<int, 3>> row_positions_unsorted;
    SimpleArrayBuffer<std::array<int, 3>> row_positions_sorted;
    SimpleArrayBuffer<RowBlock<12 * 12>> merged_rows;
};

// Accumulate species without clearing existing contributions. Particles must
// already be grouped by cell. Neutral species do not contribute.
// L = (dt/2) dJ_mid/dE_mid: the same Linear Yee weights and magnetic rotation
// as predict_current/predict_velocity are required for the discrete work balance.
// Select floor-based support and evaluate weights in logical cell coordinates;
// translate to ghost storage afterwards, so the halo does not affect weights.
template <ShapeType ST = SHAPE>
void accumulate_species(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, double dt,
                        RowBlocks& row_blocks);

// Replace the matrix from accumulated blocks; boundary constraints belong to
// the scheme and are applied afterwards to L, the RHS and the solved field.
void assemble_matrix(Operator& matrix, const RowBlocks& row_blocks, Workspace& workspace);
}   // namespace ecsim_response
