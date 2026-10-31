#pragma once

#include "boundary_conditions.h"
#include "grid/domain.h"
#include "thread_partitioned_matrix.h"

// Yee placement: E=(x+1/2,y,z), (x,y+1/2,z), (x,y,z+1/2);
// B=(x,y+1/2,z+1/2), (x+1/2,y,z+1/2), (x+1/2,y+1/2,z).
// Ordinary sparse stencils for constructing a scheme's field equation.
// Temporary during initialization; no dt or particle response is included.
struct MaxwellStencils {
    Operator curlE;
    Operator curlB;
    Operator divE;
};

MaxwellStencils assemble_maxwell_stencils(const Domain& domain, BoundaryConditionHandler& boundaries);

// Spatial Yee operators, independent of dt and particle response. Rebuild after
// changing the grid or boundary geometry; reconstructed rather than checkpointed.
struct MaxwellOperators {
    void init(MaxwellStencils&& stencils);

    ThreadPartitionedSparseMatrix<double> curlE;   // electric edges -> magnetic faces
    ThreadPartitionedSparseMatrix<double> curlB;   // magnetic faces -> electric edges
    Operator divE;                                 // electric/current grid -> charge nodes

    // Faraday: B^{n+1} = B^n - dt/2 curl(E^n + E^{n+1}).
    // The caller applies the magnetic boundary conditions to the result.
    void advance_magnetic_field(Field3d& Bn, const Field3d& B, const Field3d& E, const Field3d& En, double dt) const;
};
