#pragma once

#include "core/vector3.h"
#include "coupling/shape.h"
#include "coupling/stencil.h"
#include "grid/grid_node_view.h"
#include "include/config.h"
#include "utils/containers.h"
#include "utils/util.h"
// Logical-node gather. Two one-dimensional weight sets cover all Yee components.
// Independent raw reference kernels live in tests/unit/raw_gather_reference.h.
template <int N>
Vector3R interpolateE(const GridNodeView<const Field3d>& fieldE, const Stencil3D<N>& no,
                      const Stencil3D<N>& sh) {
    constexpr Axis X = Axis::X;
    constexpr Axis Y = Axis::Y;
    constexpr Axis Z = Axis::Z;
    double x_sum = 0, y_sum = 0, z_sum = 0;

#pragma omp simd collapse(3) reduction(+ : x_sum, y_sum, z_sum)
    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                const int ix = no.node(x, X);
                const int iy = no.node(y, Y);
                const int iz = no.node(z, Z);
                const int ix05 = sh.node(x, X);
                const int iy05 = sh.node(y, Y);
                const int iz05 = sh.node(z, Z);
                x_sum += fieldE(ix05, iy, iz, X) * sh(x, X) * no(y, Y) * no(z, Z);
                y_sum += fieldE(ix, iy05, iz, Y) * no(x, X) * sh(y, Y) * no(z, Z);
                z_sum += fieldE(ix, iy, iz05, Z) * no(x, X) * no(y, Y) * sh(z, Z);
            }
        }
    }
    return Vector3R(x_sum, y_sum, z_sum);
}

template <int N>
Vector3R interpolateB(const GridNodeView<const Field3d>& fieldB, const Stencil3D<N>& no,
                      const Stencil3D<N>& sh) {
    constexpr Axis X = Axis::X;
    constexpr Axis Y = Axis::Y;
    constexpr Axis Z = Axis::Z;
    double x_sum = 0, y_sum = 0, z_sum = 0;

#pragma omp simd collapse(3) reduction(+ : x_sum, y_sum, z_sum)
    for (int z = 0; z < N; ++z) {
        for (int y = 0; y < N; ++y) {
            for (int x = 0; x < N; ++x) {
                const int ix = no.node(x, X);
                const int iy = no.node(y, Y);
                const int iz = no.node(z, Z);
                const int ix05 = sh.node(x, X);
                const int iy05 = sh.node(y, Y);
                const int iz05 = sh.node(z, Z);
                x_sum += fieldB(ix, iy05, iz05, X) * no(x, X) * sh(y, Y) * sh(z, Z);
                y_sum += fieldB(ix05, iy, iz05, Y) * sh(x, X) * no(y, Y) * sh(z, Z);
                z_sum += fieldB(ix05, iy05, iz, Z) * sh(x, X) * sh(y, Y) * no(z, Z);
            }
        }
    }
    return Vector3R(x_sum, y_sum, z_sum);
}

// Factored Linear kernel: preserve the eight-corner evaluation order of the
// optimized raw gather while taking origins and weights from logical stencils.
template <int Component>
inline double gather_linear_component(const GridNodeView<const Field3d>& field, const Stencil3D<2>& sx,
                                      const Stencil3D<2>& sy, const Stencil3D<2>& sz) {
    const int i = sx.origin.x(), j = sy.origin.y(), k = sz.origin.z();
    const double x0 = sx(0, Axis::X), x1 = sx(1, Axis::X);
    const double y0 = sy(0, Axis::Y), y1 = sy(1, Axis::Y);
    const double z0 = sz(0, Axis::Z), z1 = sz(1, Axis::Z);
    return x0 * (y0 * (z0 * field(i, j, k, Component) + z1 * field(i, j, k + 1, Component)) +
                 y1 * (z0 * field(i, j + 1, k, Component) + z1 * field(i, j + 1, k + 1, Component))) +
           x1 * (y0 * (z0 * field(i + 1, j, k, Component) + z1 * field(i + 1, j, k + 1, Component)) +
                 y1 * (z0 * field(i + 1, j + 1, k, Component) + z1 * field(i + 1, j + 1, k + 1, Component)));
}
inline Vector3R interpolateE(const GridNodeView<const Field3d>& field, const Stencil3D<2>& no,
                             const Stencil3D<2>& sh) {
    return {gather_linear_component<0>(field, sh, no, no), gather_linear_component<1>(field, no, sh, no),
            gather_linear_component<2>(field, no, no, sh)};
}
inline Vector3R interpolateB(const GridNodeView<const Field3d>& field, const Stencil3D<2>& no,
                             const Stencil3D<2>& sh) {
    return {gather_linear_component<0>(field, no, sh, sh), gather_linear_component<1>(field, sh, no, sh),
            gather_linear_component<2>(field, sh, sh, no)};
}
