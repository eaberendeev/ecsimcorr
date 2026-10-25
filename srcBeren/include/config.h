#pragma once

#include "types.h"

#define USE_ECSIM_CORRECTION true

// ECSIM response (gather/predicted current/mass matrix) must use the same shape.
// Corrected charge/current deposition is a separate, quadratic policy.
inline constexpr ShapeType SHAPE = ShapeType::Linear;
inline constexpr ShapeType SHAPE_CH = ShapeType::Quadratic;

#define SHAPE_SIZE  2
#define GHOST_CELLS 1
#define GHOST_NODES 3

#define LMAT_MAX_ELEMENTS_PER_ROW 130
#define LMAT_VALUE_TOLERANCE      1.e-16

#define SLE_SOLVER_MAX_ITERATIONS 300
#define SLE_SOLVER_TOLERANCE      1.e-10
