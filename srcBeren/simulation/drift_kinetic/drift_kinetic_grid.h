#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "coupling/shape.h"
#include "coupling/stencil.h"
#include "drift_kinetic_push.h"
#include "grid/domain.h"
#include "utils/containers.h"

namespace driftKinetic {

using VectorField = Array3D<Vector3R>;

// Cell-centred electromagnetic fields used by the drift-kinetic model.
// Fields live on storage cells; interpolation and deposition use a linear
// two-node kernel. Gradients and the curl use centred differences.
class DriftKineticGrid {
   public:
    explicit DriftKineticGrid(const Domain& domain, std::array<bool, 3> periodic = {false, false, false})
        : domain_(domain), periodic_(periodic), electric_(domain.size()), magnetic_(domain.size()),
          gradB_(domain.size()), rotB_(domain.size()) {
    }

    VectorField& electric() {
        return electric_;
    }
    const VectorField& electric() const {
        return electric_;
    }
    VectorField& magnetic() {
        derivatives_valid_ = false;
        return magnetic_;
    }
    const VectorField& magnetic() const {
        return magnetic_;
    }

    Fields sample(const Vector3R& position) const {
        if (!derivatives_valid_)
            update_derivatives();
        const Vector3R centred = domain_.to_cell_coordinates(position) - Vector3R(0.5);
        const auto shape = make_point_stencil<ShapeType::Quadratic>(centred);

        Vector3R E(0, 0, 0), B(0, 0, 0), gradB(0, 0, 0), rotB(0, 0, 0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) {
                    const Vector3I cell{shape.origin.x() + i + GHOST_CELLS, shape.origin.y() + j + GHOST_CELLS,
                                        shape.origin.z() + k + GHOST_CELLS};
                    const double weight = shape(i, Axis::X) * shape(j, Axis::Y) * shape(k, Axis::Z);
                    E += weight * at(electric_, cell);
                    B += weight * at(magnetic_, cell);
                    gradB += weight * at(gradB_, cell);
                    rotB += weight * at(rotB_, cell);
                }

        const double lenB = B.length();
        if (!(lenB > 0.0) || !std::isfinite(lenB))
            throw std::runtime_error("drift-kinetic grid has a zero or non-finite magnetic field");
        return {E, lenB, B / lenB, gradB, rotB};
    }

    void deposit_current(VectorField& current, const Vector3R& position, const Vector3R& velocity,
                         double weight) const {
        const Vector3I cell = domain_.storage_cell_index(position);
        deposit(current, cell, position, velocity * weight);
    }

    void deposit_magnetization(VectorField& magnetization, const Vector3R& position, const Vector3R& b, double mu,
                               double weight) const {
        const Vector3I cell = domain_.storage_cell_index(position);
        deposit(magnetization, cell, position, -mu * weight * b);
    }

    // Deposits the current along a straight guiding-center segment. The
    // segment is split at cell boundaries; the longitudinal current fraction
    // of each piece is proportional to its parameter length. Magnetization is
    // deposited half at each segment endpoint, matching the midpoint rule.
    void deposit_segment(VectorField& current, VectorField& magnetization, const Vector3R& start,
                         const Vector3R& end, const Vector3R& velocity, const Vector3R& b, double charge, double mu,
                         double weight) const {
        const Vector3R a = domain_.to_cell_coordinates(start), z = domain_.to_cell_coordinates(end);
        std::vector<double> cuts{0.0, 1.0};
        for (int axis = 0; axis < 3; ++axis) {
            const double lo = std::min(a[axis], z[axis]), hi = std::max(a[axis], z[axis]);
            for (double plane = std::floor(lo) + 1.0; plane < hi; plane += 1.0) {
                const double t = (plane - a[axis]) / (z[axis] - a[axis]);
                if (t > 0.0 && t < 1.0)
                    cuts.push_back(t);
            }
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        for (size_t i = 1; i < cuts.size(); ++i) {
            const double t0 = cuts[i - 1], t1 = cuts[i], fraction = t1 - t0;
            if (!(fraction > 0.0))
                continue;
            const Vector3R r0 = start + t0 * (end - start), r1 = start + t1 * (end - start);
            const Vector3R midpoint = 0.5 * (r0 + r1);
            deposit_current(current, midpoint, velocity * (charge * fraction), weight);
            deposit_magnetization(magnetization, r0, b, mu, 0.5 * weight * fraction);
            deposit_magnetization(magnetization, r1, b, mu, 0.5 * weight * fraction);
        }
    }

    // Recompute the cached |B| gradient and curl B once per field update.
    // sample() only interpolates these arrays, not the finite differences.
    void update_derivatives() const {
        for (int i = 0; i < domain_.size().x(); ++i)
            for (int j = 0; j < domain_.size().y(); ++j)
                for (int k = 0; k < domain_.size().z(); ++k) {
                    gradB_(i, j, k) = gradient(Vector3I(i, j, k));
                    rotB_(i, j, k) = curl(Vector3I(i, j, k));
                }
        derivatives_valid_ = true;
    }

    // Cell-centred curl on the same periodic/clamped stencil as sample().
    // Used for the magnetization current curl(M) in the Maxwell source.
    VectorField curl(const VectorField& field) const {
        VectorField result(domain_.size());
        for (int i = 0; i < domain_.size().x(); ++i)
            for (int j = 0; j < domain_.size().y(); ++j)
                for (int k = 0; k < domain_.size().z(); ++k)
                    result(i, j, k) = curl(field, Vector3I(i, j, k));
        return result;
    }

   private:
    int wrap(int index, int axis) const {
        const int ghost = GHOST_CELLS, n = domain_.num_cells(axis);
        if (!periodic_[axis])
            return std::clamp(index, 0, domain_.grid.size(axis) - 1);
        return ((index - ghost) % n + n) % n + ghost;
    }

    Vector3R at(const VectorField& field, const Vector3I& cell) const {
        return field(wrap(cell.x(), 0), wrap(cell.y(), 1), wrap(cell.z(), 2));
    }

    Vector3R gradient(const Vector3I& cell) const {
        Vector3R result(0, 0, 0);
        for (int axis = 0; axis < 3; ++axis) {
            Vector3I lower = cell, upper = cell;
            --lower[axis];
            ++upper[axis];
            const double step = 2.0 * domain_.cell_size(axis);
            result[axis] = (at(magnetic_, upper).length() - at(magnetic_, lower).length()) / step;
        }
        return result;
    }

    Vector3R curl(const VectorField& field, const Vector3I& cell) const {
        auto derivative = [&](int component, int axis) {
            Vector3I lower = cell, upper = cell;
            --lower[axis];
            ++upper[axis];
            return (at(field, upper)[component] - at(field, lower)[component]) / (2.0 * domain_.cell_size(axis));
        };
        return {derivative(2, 1) - derivative(1, 2),
                derivative(0, 2) - derivative(2, 0),
                derivative(1, 0) - derivative(0, 1)};
    }

    Vector3R curl(const Vector3I& cell) const {
        return curl(magnetic_, cell);
    }

    void deposit(VectorField& field, const Vector3I& cell, const Vector3R& position, const Vector3R& value) const {
        const Vector3R centred = domain_.to_cell_coordinates(position) - Vector3R(0.5);
        const auto shape = make_point_stencil<ShapeType::Quadratic>(centred);
        (void) cell;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) {
                    const Vector3I target{shape.origin.x() + i + GHOST_CELLS, shape.origin.y() + j + GHOST_CELLS,
                                          shape.origin.z() + k + GHOST_CELLS};
                    const double weight = shape(i, Axis::X) * shape(j, Axis::Y) * shape(k, Axis::Z);
                    field(wrap(target.x(), 0), wrap(target.y(), 1), wrap(target.z(), 2)) += weight * value;
                }
    }

    const Domain& domain_;
    std::array<bool, 3> periodic_;
    VectorField electric_;
    VectorField magnetic_;
    mutable VectorField gradB_;
    mutable VectorField rotB_;
    mutable bool derivatives_valid_ = false;
};

}   // namespace driftKinetic
