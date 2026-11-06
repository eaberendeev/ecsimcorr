#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "coupling/shape.h"
#include "utils/containers.h"
#include "grid/domain.h"

// Lowest-order Whitney forms on the Yee mesh, as in the symplectic scheme.
// Charge: S1(x)S1(y)S1(z). E_x/J_x: S0(x)S1(y)S1(z).
// B_x: S1(x)S0(y)S0(z); other components follow by cyclic permutation.
// These are point weights. Exact continuity additionally needs a path integral.
namespace ecsimWhitney {
inline int periodic_index(int storage_index, int axis, const Domain& domain) {
    const int n = domain.num_cells(axis), g = domain.grid.ghost_cells();
    const int i = storage_index - g;
    return (i % n + n) % n + g;
}

// Periodic axes deposit on their owner node; non-periodic axes keep the
// physical/ghost node so the wall layer can be constrained normally. Clamping
// is a safety net for a point exactly on (or just outside) a wall where the
// upper tensor-product node would otherwise be one past the ghost layer.
inline int shape_index(int storage_index, int axis, const Domain& domain) {
    if (!domain.is_periodic(axis))
        return std::clamp(storage_index, 0, domain.grid.size(axis) - 1);
    return periodic_index(storage_index, axis, domain);
}

// Cell and local fraction of a position. At a non-periodic high wall u=n,
// choosing cell n-1 with fraction 1 keeps both nodes inside the stored grid
// and places the particle exactly on the wall edge.
inline void shape_cell_fraction(const Vector3R& position, const Domain& domain, Vector3I& cell, Vector3R& f) {
    const auto u = domain.to_cell_coordinates(position);
    for (int a = 0; a < 3; ++a) {
        if (!std::isfinite(u[a]) || std::abs(u[a]) > 1e8)
            throw std::runtime_error("non-finite or unbounded Whitney position");
        cell[a] = static_cast<int>(std::floor(u[a]));
        f[a] = u[a] - cell[a];
        if (!domain.is_periodic(a)) {
            const int n = domain.num_cells(a);
            if (cell[a] < 0) {
                cell[a] = 0;
                f[a] = 0;
            } else if (cell[a] >= n) {
                cell[a] = n - 1;
                f[a] = 1;
            }
        }
    }
}

// Copy aliases; never sum an already folded deposit a second time. On a
// non-periodic axis shape_index is the identity, so this only folds periodic
// seams.
inline void sync_periodic(Field3d& field, const Domain& domain) {
    const auto size = domain.size();
    for (int i = 0; i < size[0]; ++i)
        for (int j = 0; j < size[1]; ++j)
            for (int k = 0; k < size[2]; ++k)
                for (int c = 0; c < field.nd(); ++c)
                    field(i, j, k, c) =
                        field(shape_index(i, 0, domain), shape_index(j, 1, domain), shape_index(k, 2, domain), c);
}

struct EdgeWeight {
    int index;
    int component;
    double weight;
};

inline std::array<EdgeWeight, 12> electric_weights(const Vector3R& position, const Domain& domain) {
    Vector3I cell;
    Vector3R f;
    shape_cell_fraction(position, domain, cell, f);
    std::array<EdgeWeight, 12> result;
    for (int c = 0; c < 3; ++c) {
        const int a = (c + 1) % 3, b = (c + 2) % 3;
        for (int ia = 0; ia < 2; ++ia)
            for (int ib = 0; ib < 2; ++ib) {
                auto node = cell;
                node[a] += ia;
                node[b] += ib;
                for (int d = 0; d < 3; ++d) node[d] = shape_index(node[d] + domain.grid.ghost_cells(), d, domain);
                result[4 * c + 2 * ia + ib] = {domain.vind(node[0], node[1], node[2], c), c,
                                               (ia ? f[a] : 1 - f[a]) * (ib ? f[b] : 1 - f[b])};
            }
    }
    return result;
}

inline Vector3R gather_electric(const Field3d& field, const Vector3R& position, const Domain& domain) {
    Vector3R result(0, 0, 0);
    for (const auto& w : electric_weights(position, domain)) result[w.component] += w.weight * field[w.index];
    return result;
}

inline Vector3R gather_magnetic(const Field3d& field, const Vector3R& position, const Domain& domain) {
    Vector3I cell;
    Vector3R f, result(0, 0, 0);
    shape_cell_fraction(position, domain, cell, f);
    for (int c = 0; c < 3; ++c)
        for (int upper = 0; upper < 2; ++upper) {
            auto node = cell;
            node[c] += upper;
            for (int a = 0; a < 3; ++a) node[a] = shape_index(node[a] + domain.grid.ghost_cells(), a, domain);
            result[c] += (upper ? f[c] : 1 - f[c]) * field(node[0], node[1], node[2], c);
        }
    return result;
}

// v_bar = alpha(B) [v^n + (q dt / 2m) E_bar].
inline Vector3R magnetic_response(const Vector3R& value, const Vector3R& b) {
    return (value + value.cross(b) + b * value.dot(b)) / (1 + b.squaredNorm());
}
}   // namespace ecsimWhitney
