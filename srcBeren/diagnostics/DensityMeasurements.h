#pragma once

#include "Diagnostic.h"
#include "boundary_conditions.h"

// Charged species deposit charge; neutral species deposit particle number.
// The caller owns the result; no diagnostic grid is kept by ParticlesArray.
Field3d measure_density(const ParticlesArray& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                        ShapeType shape);
Field3d measure_charge_density(const Species& species, const Domain& domain, const BoundaryConditionHandler& boundaries,
                               ShapeType shape);
ChargeContinuityResidual measure_charge_continuity(const Field3d& before, const Field3d& after, const Field3d& current,
                                                   const Operator& divergence, double dt);
