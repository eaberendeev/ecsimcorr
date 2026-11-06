#pragma once

#include "coupling/stencil.h"
#include "grid/grid.h"

namespace ecsim_response {
// Logical Linear support for response gather and L. Keep the response kernel's
// storage-coordinate arithmetic: w1 = (coord_storage - shift) - cell_storage,
// w0 = 1 - w1. Evaluating Shape(node - coord) would round w1 differently.
// Integer truncation retains the raw kernels' support selection; particles in
// physical cells with the required lower ghost layer have nonnegative u.
inline Stencil3D<2> make_response_stencil(const Vector3R& storage_coord, const Grid& grid, double shift = 0.0) {
    Stencil3D<2> stencil;
    Vector3I storage_origin;
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        const double u = storage_coord[axis] - shift;
        storage_origin[axis] = static_cast<int>(u);
        const double fraction = u - storage_origin[axis];
        stencil.weights[static_cast<int>(axis) * 2 + 1] = fraction;
        stencil.weights[static_cast<int>(axis) * 2] = 1 - fraction;
    }
    stencil.origin = grid.to_node_index(storage_origin);
    return stencil;
}
}   // namespace ecsim_response
