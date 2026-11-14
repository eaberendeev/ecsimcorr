#pragma once

#include <array>
#include <vector>

#include "fields/maxwell_operators.h"

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
    // Sum of species responses on the unknowns, on the fixed pattern of one particle's
    // couplings united with I-M; cleared before each predictor assembly.
    // All matrices are reconstructed on restart, never checkpointed.
    Operator L;

    void clear_response();
    // I-M+L; rows without a plasma response keep the compact I-M row.
    const Operator& predictor_matrix();

   private:
    Operator fixed_columns_;              // I-M of unknown rows in fixed stored components
    std::vector<double> im_on_pattern_;   // I-M values at the entries of L
    Operator predictor_;                  // scratch of predictor_matrix()
};
