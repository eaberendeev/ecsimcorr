#pragma once

enum class Axis : int {
    X = 0,
    Y = 1,
    Z = 2,
    Comp = 3,   // Component coordinate in the (x, y, z, component) field layout.
};
// Compile-time default halo for node views and default domain construction.
// A concrete Grid still owns its actual runtime halo width.
inline constexpr int GHOST_CELLS = 1;

enum class ShapeType { NGP, Linear, Quadratic };
