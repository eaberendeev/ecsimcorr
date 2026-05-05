#pragma once

#define USE_ECSIM_CORRECTION true

#define SHAPE    ShapeType::Linear
#define SHAPE_CH ShapeType::Quadratic

#define SHAPE_SIZE  2
#define GHOST_CELLS 1
#define GHOST_NODES 3

#define LMAT_MAX_ELEMENTS_PER_ROW 130
#define LMAT_VALUE_TOLERANCE      1.e-16

#define SLE_SOLVER_MAX_ITERATIONS 300
#define SLE_SOLVER_TOLERANCE      1.e-10
