#pragma once

enum class Axis : int {
    X = 0,
    Y = 1,
    Z = 2,
    Comp = 3,   // Component coordinate in the (x, y, z, component) field layout.
};
enum class ShapeType { NGP, Linear, Quadratic };
