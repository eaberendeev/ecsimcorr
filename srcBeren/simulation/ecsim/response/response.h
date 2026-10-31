#pragma once

#include "bmatrix.h"
#include "config.h"
#include "grid/domain.h"
#include "memory.h"
#include "row_block.h"

class ParticlesArray;

namespace ecsim_response {
// Per-thread response contributions are assembly results, passed explicitly
// from species accumulation to CSR conversion. Reset once before all species.
using RowBlocks = SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>>;
void reset_rows(RowBlocks& row_blocks);

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
template <ShapeType ST = SHAPE>
void accumulate_species(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, double dt,
                        RowBlocks& row_blocks);
template <ShapeType ST = SHAPE>
void accumulate_species_reference(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, double dt,
                                  BlockMatrix& reference_blocks);

// Replace the matrix from accumulated blocks; boundary constraints belong to
// the scheme and are applied afterwards to L, the RHS and the solved field.
void assemble_matrix(Operator& matrix, const RowBlocks& row_blocks, Workspace& workspace);
void assemble_reference(Operator& matrix, const Domain& domain, const BlockMatrix& reference_blocks);
}   // namespace ecsim_response
