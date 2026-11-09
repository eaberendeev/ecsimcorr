#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "coupling/shape.h"
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
        const Vector3R u = domain_.to_cell_coordinates(position), centred = u - Vector3R(0.5);
        Vector3I base;
        Vector3R fraction;
        for (int axis = 0; axis < 3; ++axis) {
            const int raw = static_cast<int>(std::floor(centred[axis]));
            base[axis] = raw + GHOST_CELLS;
            fraction[axis] = centred[axis] - raw;
        }

        Vector3R E(0, 0, 0), B(0, 0, 0), gradB(0, 0, 0), rotB(0, 0, 0);
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k) {
                    const Vector3I cell{base.x() + i, base.y() + j, base.z() + k};
                    const double weight = (i ? fraction.x() : 1 - fraction.x()) *
                                          (j ? fraction.y() : 1 - fraction.y()) *
                                          (k ? fraction.z() : 1 - fraction.z());
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
        const Vector3R u = domain_.to_cell_coordinates(position), centred = u - Vector3R(0.5);
        Vector3I base;
        Vector3R fraction;
        for (int axis = 0; axis < 3; ++axis) {
            const int raw = static_cast<int>(std::floor(centred[axis]));
            base[axis] = raw + GHOST_CELLS;
            fraction[axis] = centred[axis] - raw;
        }
        (void) cell;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k) {
                    const Vector3I target{base.x() + i, base.y() + j, base.z() + k};
                    const double weight = (i ? fraction.x() : 1 - fraction.x()) *
                                          (j ? fraction.y() : 1 - fraction.y()) *
                                          (k ? fraction.z() : 1 - fraction.z());
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
