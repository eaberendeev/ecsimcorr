#pragma once

#include "maxwell_operators.h"

// Shared Maxwell part of the ECSIM predictor and Corr endpoint equation.
// Owned by the scheme; rebuild when dt, grid or field boundaries change.
struct EcsimFieldSystem {
    void init(const MaxwellStencils& stencils, const Domain& domain, double dt,
              const BoundaryConditionHandler& boundaries);

    Operator M;    // -dt^2/4 curlB curlE, before electric Dirichlet constraints
    Operator IM;   // I-M, with prescribed electric rows pinned to identity
    // Sum of species responses; overwritten at each predictor assembly.
    // All matrices are reconstructed on restart, never checkpointed.
    Operator L;
};
