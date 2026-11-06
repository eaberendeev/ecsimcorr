#pragma once

#include "core/indexing.h"
#include "coupling/shape.h"

// Tensor-product point weights on logical nodes; independent of ghost storage.
template <int N>
struct Stencil3D {
    static_assert(N > 0);
    alignas(64) std::array<double, 3 * N> weights;
    Vector3I origin;

    int node(int local, Axis axis) const {
        assert(local >= 0 && local < N);
        assert(static_cast<int>(axis) >= 0 && static_cast<int>(axis) < 3);
        return origin[axis] + local;
    }
    double& operator()(int local, Axis axis) {
        assert(local >= 0 && local < N);
        assert(static_cast<int>(axis) >= 0 && static_cast<int>(axis) < 3);
        return weights[static_cast<int>(axis) * N + local];
    }
    const double& operator()(int local, Axis axis) const {
        assert(local >= 0 && local < N);
        assert(static_cast<int>(axis) >= 0 && static_cast<int>(axis) < 3);
        return weights[static_cast<int>(axis) * N + local];
    }
};

template <ShapeFunction ShapeFn, int N>
Stencil3D<N> make_stencil_at(const Vector3R& coord, const Vector3I& origin, const Vector3R& shift = Vector3R(0.0)) {
    for ([[maybe_unused]] Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        assert(origin[axis] <= std::numeric_limits<int>::max() - (N - 1));
        assert(std::isfinite(coord[axis]) && std::isfinite(shift[axis]) && std::isfinite(coord[axis] - shift[axis]));
    }
    Stencil3D<N> stencil;
    stencil.origin = origin;
    const Vector3R u = coord - shift;
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z})
        for (int local = 0; local < N; ++local)
            stencil(local, axis) = ShapeFn(double(stencil.node(local, axis)) - u[axis]);
    return stencil;
}

template <ShapeType ST>
Stencil3D<ST == ShapeType::Linear ? 2 : 3> make_point_stencil(const Vector3R& coord,
                                                              const Vector3R& shift = Vector3R(0.0)) {
    static_assert(ST == ShapeType::Linear || ST == ShapeType::Quadratic);
    constexpr int N = ST == ShapeType::Linear ? 2 : 3;
    const Vector3R u = coord - shift;
    Vector3I origin;
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z})
        origin[axis] =
            checked_floor_index(u[axis] + (ST == ShapeType::Linear ? 0.0 : 0.5), ST == ShapeType::Linear ? 0 : -1);
    return make_stencil_at<ShapeKernel<ST>::weight, N>(coord, origin, shift);
}

// Four-node deposition envelope anchored to the particle's logical cell.
inline Vector3I deposition_origin(const Vector3R& coord) {
    return checked_floor_index(coord, -1);
}
inline Vector3I deposition_origin(const Vector3I& cell) {
    for (int d = 0; d < 3; ++d) assert(cell[d] > std::numeric_limits<int>::min());
    return cell - Vector3I(1);
}
template <ShapeFunction ShapeFn, int N = 4>
Stencil3D<N> make_deposition_stencil(const Vector3R& coord) {
    static_assert(N == 4);
    return make_stencil_at<ShapeFn, N>(coord, deposition_origin(coord));
}
template <ShapeFunction ShapeFn, int N>
bool stencil_support_fits(const Stencil3D<N>& stencil, const Vector3R& point) {
    static_assert(ShapeFn == Shape || ShapeFn == Shape2);
    constexpr double radius = ShapeFn == Shape2 ? 1.5 : 1.0;
    for (int d = 0; d < 3; ++d) {
        const double lower = double(stencil.origin[d]) + radius - 1.0;
        const double upper = double(stencil.origin[d]) + N - radius;
        if (!std::isfinite(point[d]) || point[d] < lower || point[d] > upper)
            return false;
    }
    return true;
}
template <ShapeFunction ShapeFn>
Vector3I trajectory_stencil_origin(const Vector3R& from, const Vector3R& to) {
    static_assert(ShapeFn == Shape || ShapeFn == Shape2);
    constexpr double shift = ShapeFn == Shape2 ? 0.5 : 0.0;
    Vector3I origin;
    for (int d = 0; d < 3; ++d) origin[d] = checked_floor_index(std::min(from[d], to[d]) + shift, -1);
    return origin;
}
