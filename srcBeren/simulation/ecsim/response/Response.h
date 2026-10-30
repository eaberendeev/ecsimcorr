#pragma once

#include "bmatrix.h"
#include "config.h"
#include "grid/domain.h"
#include "memory.h"
#include "row_block.h"

class ParticlesArray;

namespace ecsim_response {
// Scratch for the Linear ECSIM response at x^{n+1/2}, B^n + B_external.
// Simulation owns this storage for the run. Row buffers are reset before each
// assembly; reference blocks are prepared only for validation. No checkpoint
// data: the response is reconstructed from particles and fields each step.
struct Workspace {
    Workspace();
    void reset_rows();
    SimpleArrayBuffer<SimpleArrayBuffer<RowBlock<36>>> row_blocks;
    BlockMatrix reference_blocks;
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
                        Workspace& workspace);
template <ShapeType ST = SHAPE>
void accumulate_species_reference(const ParticlesArray& species, const Field3d& fieldB, const Domain& domain, double dt,
                                  Workspace& workspace);

// Replace the matrix from accumulated blocks; boundary constraints belong to
// the scheme and are applied afterwards to L, the RHS and the solved field.
void assemble_matrix(Operator& matrix, Workspace& workspace);
void assemble_reference(Operator& matrix, const Domain& domain, const Workspace& workspace);
}   // namespace ecsim_response
