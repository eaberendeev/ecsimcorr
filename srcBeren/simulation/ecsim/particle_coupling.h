// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2024, for licensing details see the LICENSE file

#pragma once

#include "particle_array.h"
#include "particles/particle_path.h"
#include "response/response.h"

// The standard ECSIM coupling keeps Jp, L and the Boris gather consistent.
// Response weights use SHAPE (Linear); Corr trajectory charge uses SHAPE_CH.
// This concrete implementation does not select Whitney or other response models.
class LinearEcsimCoupling {
   public:
    LinearEcsimCoupling();

    // Assemble L = (dt/2) dJ_mid/dE_mid from particles already grouped by cell.
    // B_full includes the prescribed magnetic field. Boundary constraints are
    // applied by the scheme after this call; species and B_full are unchanged.
    void assemble_response(Species& species, const Field3d& B_full, const Domain& domain, double dt, Operator& L);
    // Accumulate field-independent Jp without clearing existing contributions.
    void predict_current(const ParticlesArray& species, const Field3d& B_full, Field3d& J, double dt) const;
    // v^n -> v' using the same gather as the assembled response; coordinates do not move.
    void push_velocities(ParticlesArray& species, const Field3d& E_full, const Field3d& B_full, double dt) const;
    // Corr: Boris over dt, then drift/deposit over half_dt. Accumulate the
    // half-path current and predictor power without applying recorded wall fates.
    void push_and_deposit(ParticlesArray& species, const Field3d& E_full, const Field3d& B_full, double dt,
                          double half_dt, Field3d& J, double& predictor_power, const PathPolicy& policy) const;

   private:
    void prepare_reference_blocks(const Species& species, const Domain& domain);

    // Assembly scratch is overwritten on every call; retain allocations only.
    // The validation counter belongs to this coupling instance. Neither it nor
    // these buffers are physical state or part of the current recovery format.
    ecsim_response::RowBlocks response_row_blocks_;
    BlockMatrix reference_response_blocks_;
    ecsim_response::Workspace response_workspace_;
    int assembly_check_counter_ = 0;
};
