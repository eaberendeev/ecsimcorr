#pragma once

#include "core/indexing.h"
#include "coupling/shape.h"

// Tensor-product point weights on logical nodes; independent of ghost storage.
template <int N>
struct Stencil3D {
    static_assert(N > 0);
    Vector3I origin;
    std::array<double, 3 * N> weights;

    int node(int local, Axis axis) const {
        assert(local >= 0 && local < N);
        assert(static_cast<int>(axis) >= 0 && static_cast<int>(axis) < 3);
        return origin[axis] + local;
    }
    double operator()(int local, Axis axis) const {
        assert(local >= 0 && local < N);
        assert(static_cast<int>(axis) >= 0 && static_cast<int>(axis) < 3);
        return weights[static_cast<int>(axis) * N + local];
    }
};

template <auto ShapeFn, int N>
Stencil3D<N> make_stencil_at(const Vector3R& coord, const Vector3I& origin, const Vector3R& shift = Vector3R(0.0)) {
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        if (origin[axis] > std::numeric_limits<int>::max() - (N - 1))
            throw std::out_of_range("Stencil support exceeds integer range");
        if (!std::isfinite(coord[axis]) || !std::isfinite(shift[axis]) || !std::isfinite(coord[axis] - shift[axis]))
            throw std::out_of_range("Stencil coordinates must be finite");
    }
    Stencil3D<N> stencil{origin, {}};
    const Vector3R u = coord - shift;
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z})
        for (int local = 0; local < N; ++local)
            stencil.weights[static_cast<int>(axis) * N + local] = ShapeFn(double(stencil.node(local, axis)) - u[axis]);
    return stencil;
}

template <ShapeType ST>
auto make_point_stencil(const Vector3R& coord, const Vector3R& shift = Vector3R(0.0)) {
    static_assert(ST == ShapeType::Linear || ST == ShapeType::Quadratic);
    constexpr int N = ST == ShapeType::Linear ? 2 : 3;
    Stencil3D<N> stencil;
    const Vector3R u = coord - shift;
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        // The checked floor already validates the shifted coordinate. Build
        // this axis directly instead of validating it again in make_stencil_at.
        const int origin =
            checked_floor_index(u[axis] + (ST == ShapeType::Linear ? 0.0 : 0.5), ST == ShapeType::Linear ? 0 : -1);
        if (origin > std::numeric_limits<int>::max() - (N - 1))
            throw std::out_of_range("Stencil support exceeds integer range");
        stencil.origin[axis] = origin;
        for (int local = 0; local < N; ++local)
            stencil.weights[static_cast<int>(axis) * N + local] =
                ShapeKernel<ST>::weight(double(origin + local) - u[axis]);
    }
    return stencil;
}
