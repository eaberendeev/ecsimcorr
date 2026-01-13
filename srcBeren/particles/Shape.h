#pragma once

#include <array>

#include "vector3.h"

using ShapeFunction = double (*)(const double&);

template <auto ShapeFn, int SMAX>
struct ParticleShape {
    static_assert(SMAX > 0, "SMAX must be positive");

    alignas(64) std::array<double, 3 * SMAX> w_;
    Vector3I base_;
    Vector3I start_;      // start_ = base_ - ghost_cells
    Vector3R shift_;   // shift base from coord (in cell units)

    inline int idx(int i, int comp) const { return i + SMAX * comp; }

    const double& operator()(int i, int comp) const { return w_[idx(i, comp)]; }
    inline void fill_from_normalized(const Vector3R& coord,
                                     const Vector3R& shift = {0, 0,
                                                              0}) noexcept {
        int ghost_cells = SMAX / 2 - 1;
        return fill_from_normalized(coord, ghost_cells, shift);
    }
    // Заполнить по нормализованным координатам (xx = pos.x / cellSize и т.п.)
    inline void fill_from_normalized(const Vector3R& coord_, int ghost_cells,
                                     const Vector3R& shift) noexcept {
        Vector3R coord = coord_ - shift;
        base_.x() = floor(coord.x() );
        base_.y() = floor(coord.y() );
        base_.z() = floor(coord.z() );
        start_ = base_ - Vector3I(ghost_cells);
        shift_ = shift;

        for (int n = 0; n < SMAX; ++n) {
            const double argx = -coord.x() + double(start_.x() + n);
            const double argy = -coord.y() + double(start_.y() + n);
            const double argz = -coord.z() + double(start_.z() + n);
            w_[idx(n, Axis::X)] = ShapeFn(argx);
            w_[idx(n, Axis::Y)] = ShapeFn(argy);
            w_[idx(n, Axis::Z)] = ShapeFn(argz);
        }
    }
    // Заполнить по нормализованным координатам (xx = pos.x / cellSize и т.п.)
    inline void fill_from_normalized(const Vector3R& coord, const Vector3I& base,
                                     int ghost_cells, const Vector3R& shift = {0,0,0} ) noexcept {
        base_ = base;
        start_ = base_ - Vector3I(ghost_cells);
        shift_ = shift;

        for (int n = 0; n < SMAX; ++n) {
            const double argx = -coord.x() + double(start_.x() + n);
            const double argy = -coord.y() + double(start_.y() + n);
            const double argz = -coord.z() + double(start_.z() + n);
            w_[idx(n, Axis::X)] = ShapeFn(argx);
            w_[idx(n, Axis::Y)] = ShapeFn(argy);
            w_[idx(n, Axis::Z)] = ShapeFn(argz);
        }
    }
    inline void fill_zero() noexcept {
        base_ = {0, 0, 0};
        start_ = {0, 0, 0};
        shift_ = {0, 0, 0};
        std::fill(w_.begin(), w_.end(), 0.0);
    }
};

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

inline double Shape3(const double& dist) {
    double d = fabs(dist);
    if (d <= 1.)
        return (2. / 3. - d * d + 0.5 * d * d * d);
    else if (d < 2)
        return ((2. - d) * (2. - d) * (2. - d) / 6.);
    else
        return 0.;
}

inline double Shape4(const double& dist) {
    double d = fabs(dist);
    if (d <= 0.5)
        return (115. / 192 - 0.625 * d * d + 0.25 * d * d * d * d);
    else if (d <= 1.5)
        return (55. + 20. * d - 120. * d * d + 80. * d * d * d -
                16. * d * d * d * d) /
               96.;
    else if (d < 2.5)
        return (5. - 2. * d) * (5. - 2. * d) * (5. - 2. * d) * (5. - 2. * d) /
               384.;
    else
        return 0.;
}

template <auto ShapeFn>
constexpr int shape_width() {
    return 2;
}

template <>
constexpr int shape_width<Shape>() {
    return 2;
}  // 2 nodes
template <>
constexpr int shape_width<Shape2>() {
    return 2;
}   // Shape2: support up to 1.5 -> 2 nodes
template <>
constexpr int shape_width<Shape3>() {
    return 3;
}