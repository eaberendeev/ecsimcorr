#pragma once

enum class Axis : int {
    X = 0,
    Y = 1,
    Z = 2,
    Comp = 3,   // Component coordinate in the (x, y, z, component) field layout.
};
#ifndef BEREN_GHOST_CELLS
#define BEREN_GHOST_CELLS 1
#endif
// One storage halo width for the entire simulation; stencil geometry is independent.
inline constexpr int GHOST_CELLS = BEREN_GHOST_CELLS;
static_assert(GHOST_CELLS >= 1, "The simulation requires at least one ghost layer");

enum class ShapeType { NGP, Linear, Quadratic };
