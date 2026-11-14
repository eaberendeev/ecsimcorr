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
    // Across-node SIMD helps the four-node quadratic envelope. Linear weights
    // are faster with the axis-major loop (measured separately in Release).
    if constexpr (ShapeFn == Shape2 && N == 4) {
#pragma omp simd
        for (int local = 0; local < N; ++local) {
            stencil(local, Axis::X) = ShapeFn(double(stencil.node(local, Axis::X)) - u.x());
            stencil(local, Axis::Y) = ShapeFn(double(stencil.node(local, Axis::Y)) - u.y());
            stencil(local, Axis::Z) = ShapeFn(double(stencil.node(local, Axis::Z)) - u.z());
        }
    } else {
        for (Axis axis : {Axis::X, Axis::Y, Axis::Z})
            for (int local = 0; local < N; ++local)
                stencil(local, axis) = ShapeFn(double(stencil.node(local, axis)) - u[axis]);
    }
    return stencil;
}

// Select half-shifted Linear support inside an already known logical cell.
// Compare u-.5 itself: comparing (u-cell)<.5 changes rounding at negative half nodes.
inline Vector3I linear_half_shifted_origin(const Vector3R& coord, const Vector3I& cell) {
    assert(checked_floor_index(coord) == cell);
    const Vector3R shifted = coord - Vector3R(0.5);
    const Vector3I offset(shifted.x() < double(cell.x()), shifted.y() < double(cell.y()),
                          shifted.z() < double(cell.z()));
    for ([[maybe_unused]] Axis axis : {Axis::X, Axis::Y, Axis::Z})
        assert(offset[axis] == 0 || cell[axis] > std::numeric_limits<int>::min());
    return cell - offset;
}

// Linear point support is already known: avoid abs/clamping and repeated floor.
// coord-shift must lie between origin and origin+1, including both endpoints.
inline Stencil3D<2> make_linear_stencil_at(const Vector3R& coord, const Vector3I& origin,
                                           const Vector3R& shift = Vector3R(0.0)) {
    Stencil3D<2> stencil;
    stencil.origin = origin;
    const Vector3R u = coord - shift;
    for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
        const double fraction = u[axis] - origin[axis];
        assert(std::isfinite(fraction) && fraction >= 0.0 && fraction <= 1.0);
        assert(origin[axis] < std::numeric_limits<int>::max());
        stencil(0, axis) = 1.0 - fraction;
        stencil(1, axis) = fraction;
    }
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
        assert(std::isfinite(point[d]));
        if (point[d] < lower || point[d] > upper)
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
