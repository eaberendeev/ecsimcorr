#pragma once

#include <array>

#include "Vec.h"

using ShapeFunction = double (*)(const double&);

template <auto ShapeFn, int SMAX>
struct ParticleShape {
    static_assert(SMAX > 0, "SMAX must be positive");

    alignas(64) std::array<double, 3 * SMAX> w_;
    int3 base_;
    int3 start_;      // start_ = base_ - ghost_cells
    double3 shift_;   // shift base from coord (in cell units)

    inline int idx(int i, int comp) const { return i + SMAX * comp; }

    const double& operator()(int i, int comp) const { return w_[idx(i, comp)]; }

    // Заполнить по нормализованным координатам (xx = pos.x / cellSize и т.п.)
    inline void fill_from_normalized(const double3& coord, int ghost_cells,
                                     const double3& shift = {0, 0,
                                                             0}) noexcept {
        base_.x() = static_cast<int>(coord.x() + shift.x());
        base_.y() = static_cast<int>(coord.y() + shift.y());
        base_.z() = static_cast<int>(coord.z() + shift.z());
        start_ = base_ - int3(ghost_cells);
        shift_ = shift;

        for (int n = 0; n < SMAX; ++n) {
            const double argx = -coord.x() + double(start_.x() + n);
            const double argy = -coord.y() + double(start_.y() + n);
            const double argz = -coord.z() + double(start_.z() + n);
            w_[idx(n, Axis::X)] = ShapeFn(argx);
            w_[idx(n, Axis::Y)] = ShapeFn(argy);
            w_[idx(n, Axis::Z)] = ShapeFn(argz);
            std::cout << "argx " << argx << " " << ShapeFn(argx) << std::endl;
            std::cout << idx(n, Axis::X)  << " " << w_[idx(n, Axis::X)] << " " 
            << idx(n, Axis::Y)  << " " << w_[idx(n, Axis::Y)] << " "
            << idx(n, Axis::Z)  << " " << w_[idx(n, Axis::Z)] << std::endl;
        }
    }
    // Заполнить по нормализованным координатам (xx = pos.x / cellSize и т.п.)
    inline void fill_from_normalized(const double3& coord, const int3& base,
                                     int ghost_cells, const double3& shift = {0,0,0} ) noexcept {
        base_ = base;
        start_ = base_ - int3(ghost_cells);
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
