#pragma once

#include <vector>

#include "fields/maxwell_operators.h"

// Shared Maxwell part of the ECSIM predictor and Corr endpoint equation.
// Owned by the scheme; rebuild when dt, grid or field boundaries change.
// The equations are solved for one unknown per physical E component: the
// periodic images of a node share the unknown of its owner in [g, g+N).
// Fields keep the full storage: restrict() takes the owner values of a field,
// prolong() copies unknowns to every image.
struct EcsimFieldSystem {
    void init(const MaxwellStencils& stencils, const Domain& domain, double dt,
              const FieldBoundaryConditions& boundaries);

    int unknowns() const {
        return std::ssize(storage_index);
    }

    Field3d restrict(const Field3d& field) const;
    void prolong(const Field3d& unknowns, Field3d& field) const;

    std::vector<int> unknown_index;   // stored E component -> unknown
    std::vector<int> storage_index;   // unknown -> stored E component of the owner

    Operator M;    // -dt^2/4 curlB curlE on storage, for right-hand sides
    Operator IM;   // I-M on the unknowns, prescribed electric rows pinned to identity
    // Sum of species responses on the unknowns; overwritten at each predictor assembly.
    // All matrices are reconstructed on restart, never checkpointed.
    Operator L;
};
