#pragma once

#include <array>
#include <vector>

#include "fields/maxwell_operators.h"
#include "simulation/ecsim/response/response.h"
#include "utils/containers.h"
#include "utils/memory.h"

// Shared Maxwell part of the ECSIM predictor and Corr endpoint equation.
// Owned by the scheme; rebuild when dt, grid or field boundaries change.
// The equations are solved for one unknown per physical E component: the
// periodic images of a node share the unknown of its owner in [g, g+N).
// Components whose I-M row is the identity and which are not prescribed (nodes
// outside the field domain) are fixed: they keep their right-hand side value
// and carry no plasma response, so they are not unknowns.
// Fields keep the full storage: restrict() takes the unknowns of a field,
// prolong() copies them to every image.
struct EcsimFieldSystem {
    void init(const MaxwellStencils& stencils, const Domain& domain, double dt,
              const FieldBoundaryConditions& boundaries);

    int unknowns() const {
        return std::ssize(storage_index);
    }

    Field3d restrict(const Field3d& field) const;
    // Right-hand side on the unknowns: owner rows less the I-M columns of fixed components.
    Field3d restrict_rhs(const Field3d& rhs) const;
    // Fixed components take their value from `fixed` (the storage right-hand side), or zero.
    void prolong(const Field3d& unknowns, const Field3d& fixed, Field3d& field) const;
    void prolong(const Field3d& unknowns, Field3d& field) const;

    std::vector<int> unknown_index;   // stored E component -> unknown, -1 if fixed
    std::vector<int> storage_index;   // unknown -> stored E component of the owner

    Operator M;    // -dt^2/4 curlB curlE on storage, for right-hand sides
    Operator IM;   // I-M on the unknowns, prescribed electric rows pinned to identity
    // Sum of species responses on the unknowns, zeroed before each predictor
    // assembly. Its pattern is that of I-M, united with the couplings of one particle
    // in the rows near the covered cells (cover_cells()). All matrices are
    // reconstructed on restart, never checkpointed.
    Operator L;
    // Positions of the cell-block couplings in the rows of L.
    ecsim_response::ResponseLayout response_layout;

    // Give response rows to the cells holding particles (storage cell order). The
    // pattern of L is rebuilt only when an occupied cell is not covered yet; it then
    // covers the occupied cells and a margin around them, so a plasma that keeps its
    // region is never rebuilt. A rebuild soon after the previous one doubles the
    // margin (a growing plasma), a later one halves it back. Above
    // full_cover_fraction of the cells the whole box is covered and L keeps the
    // full pattern. Returns true on a rebuild.
    static constexpr int min_cover_margin = 2;
    static constexpr int max_cover_margin = 16;
    static constexpr int frequent_rebuild_covers = 8;
    static constexpr double full_cover_fraction = 0.7;
    // `occupied` has one flag per storage cell (ParticlesArray cell order).
    bool cover_cells(const Array3D<int8_t>& occupied);
    void zero_response();
    // I-M+L without the response entries no particle reached; rows without a plasma
    // response reduce to the compact I-M row. Assembled into a member matrix on
    // each call: the reference stays valid until the next call.
    const Operator& predictor_matrix();

   private:
    void build_response(const std::vector<int8_t>& active_rows);

    Domain domain_;
    std::array<bool, 3> periodic_;
    Array3D<int8_t> covered_cells_;
    Array3D<int8_t> node_cells_;       // in cover_cells(): cells whose nodes get response rows
    Array3D<int8_t> dilate_scratch_;   // in dilate(): the flags between its passes along y and z
    int cover_margin_ = min_cover_margin;
    int covers_since_rebuild_ = 0;   // cover_cells() calls since the last rebuild
    Operator fixed_columns_;         // I-M of unknown rows in fixed stored components
    // Nonzero I-M entries: positions in L and values, rows start at im_start_.
    SmartPtr<int> im_start_;
    SmartPtr<int> im_position_;
    SmartPtr<double> im_value_;
    SmartPtr<int> chunk_offsets_;   // scratch of predictor_matrix()
    Operator predictor_;            // scratch of predictor_matrix()
};
