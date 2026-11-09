#pragma once

#include <array>
#include <cassert>

#include "core/types.h"
#include "core/vector3.h"

using ShapeFunction = double (*)(const double&);

inline double Shape(const double& dist) {
    double d = fabs(dist);
    if (d < 1.)
        return (1. - d);
    else
        return 0.;
}

inline double Shape2(const double& dist) {
    double d = fabs(dist);
    if (d <= 0.5)
        return (0.75 - d * d);
    else if (d < 1.5)
        return ((d - 1.5) * (d - 1.5) * 0.5);
    else
        return 0.;
}

// Compile-time dispatch only. Both deposit kernels use a four-node envelope
// for a trajectory shorter than one cell, irrespective of the shape support.
template <ShapeType ST>
struct ShapeKernel {
    static_assert(ST == ShapeType::Linear || ST == ShapeType::Quadratic,
                  "density/current deposition supports Linear and Quadratic only");
    static constexpr ShapeFunction weight = ST == ShapeType::Linear ? Shape : Shape2;
    static constexpr int deposit_half_width = 2;
};
