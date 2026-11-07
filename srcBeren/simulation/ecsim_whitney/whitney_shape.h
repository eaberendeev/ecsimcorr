#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "coupling/shape.h"
#include "grid/domain.h"
#include "utils/containers.h"

// Whitney forms on the Yee mesh, as in the symplectic scheme.
//
// Lowest order (linear kernel):
//   charge: S1(x)S1(y)S1(z). E_x/J_x: S0(x)S1(y)S1(z).
//   B_x: S1(x)S0(y)S0(z); other components follow by cyclic permutation.
//
// Quadratic kernel (second order in every direction):
//   charge: S2(x)S2(y)S2(z). E_x/J_x: S1(x)S2(y)S2(z).
//   B_x: S2(x)S1(y)S1(z).
// S2 is the quadratic B-spline from coupling/shape.h (support 3 nodes).
// These are point weights. Exact continuity additionally needs a path integral.
namespace ecsimWhitney {
enum class KernelOrder { Linear, Quadratic };

constexpr int max_electric_weights = 54;   // 3 components * 2 * 3 * 3
constexpr int max_magnetic_weights = 36;   // 3 components * 3 * 2 * 2

struct EdgeWeight {
    int index;
    int component;
    double weight;
};

struct ElectricStencil {
    std::array<EdgeWeight, max_electric_weights> weight{};
    int size = 0;
};

struct MagneticStencil {
    std::array<EdgeWeight, max_magnetic_weights> weight{};
    int size = 0;
};

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

// One axis of a tensor-product shape. S0 is a single node, S1 a two-node
// hat, S2 the three-node quadratic B-spline.
struct AxisWeights {
    int base = 0;
    int count = 1;
    double weight[3] = {1, 0, 0};
};

inline AxisWeights linear_axis(double u) {
    const int cell = static_cast<int>(std::floor(u));
    const double f = u - cell;
    return {cell, 2, {1 - f, f, 0}};
}

inline AxisWeights quadratic_axis(double u) {
    const int base = static_cast<int>(std::floor(u + 0.5)) - 1;
    return {base, 3, {Shape2(u - base), Shape2(u - (base + 1)), Shape2(u - (base + 2))}};
}

inline int stencil_index(int i, int j, int k, int component, const Domain& domain) {
    const int g = domain.grid.ghost_cells();
    return domain.vind(shape_index(i + g, 0, domain), shape_index(j + g, 1, domain), shape_index(k + g, 2, domain),
                       component);
}

inline void append_weight(ElectricStencil& stencil, int i, int j, int k, int component, double weight,
                          const Domain& domain) {
    if (weight == 0)
        return;
    stencil.weight[stencil.size++] = {stencil_index(i, j, k, component, domain), component, weight};
}

inline void append_weight(MagneticStencil& stencil, int i, int j, int k, int component, double weight,
                          const Domain& domain) {
    if (weight == 0)
        return;
    stencil.weight[stencil.size++] = {stencil_index(i, j, k, component, domain), component, weight};
}

inline ElectricStencil electric_stencil_linear(const Vector3R& position, const Domain& domain) {
    Vector3I cell;
    Vector3R f;
    shape_cell_fraction(position, domain, cell, f);
    ElectricStencil result;
    for (int c = 0; c < 3; ++c) {
        const int a = (c + 1) % 3, b = (c + 2) % 3;
        for (int ia = 0; ia < 2; ++ia)
            for (int ib = 0; ib < 2; ++ib) {
                auto node = cell;
                node[a] += ia;
                node[b] += ib;
                append_weight(result, node[0], node[1], node[2], c, (ia ? f[a] : 1 - f[a]) * (ib ? f[b] : 1 - f[b]),
                              domain);
            }
    }
    return result;
}

// E_x/J_x = S1(x) S2(y) S2(z): linear along the component, quadratic across it.
inline ElectricStencil electric_stencil_quadratic(const Vector3R& position, const Domain& domain) {
    const auto u = domain.to_cell_coordinates(position);
    for (int a = 0; a < 3; ++a)
        if (!std::isfinite(u[a]) || std::abs(u[a]) > 1e8)
            throw std::runtime_error("non-finite or unbounded Whitney position");
    AxisWeights along[3], transverse[3];
    for (int a = 0; a < 3; ++a) {
        // The edge coordinate is staggered by half a cell along its own axis.
        along[a] = linear_axis(u[a] - 0.5);
        transverse[a] = quadratic_axis(u[a]);
    }
    ElectricStencil result;
    for (int c = 0; c < 3; ++c) {
        const int a = (c + 1) % 3, b = (c + 2) % 3;
        for (int ic = 0; ic < along[c].count; ++ic)
            for (int ia = 0; ia < transverse[a].count; ++ia)
                for (int ib = 0; ib < transverse[b].count; ++ib) {
                    int node[3];
                    node[c] = along[c].base + ic;
                    node[a] = transverse[a].base + ia;
                    node[b] = transverse[b].base + ib;
                    append_weight(result, node[0], node[1], node[2], c,
                                  along[c].weight[ic] * transverse[a].weight[ia] * transverse[b].weight[ib], domain);
                }
    }
    return result;
}

inline MagneticStencil magnetic_stencil_linear(const Vector3R& position, const Domain& domain) {
    Vector3I cell;
    Vector3R f;
    shape_cell_fraction(position, domain, cell, f);
    MagneticStencil result;
    for (int c = 0; c < 3; ++c)
        for (int upper = 0; upper < 2; ++upper) {
            auto node = cell;
            node[c] += upper;
            append_weight(result, node[0], node[1], node[2], c, upper ? f[c] : 1 - f[c], domain);
        }
    return result;
}

// B_x = S2(x) S1(y) S1(z): quadratic along the component, linear across it.
inline MagneticStencil magnetic_stencil_quadratic(const Vector3R& position, const Domain& domain) {
    const auto u = domain.to_cell_coordinates(position);
    for (int a = 0; a < 3; ++a)
        if (!std::isfinite(u[a]) || std::abs(u[a]) > 1e8)
            throw std::runtime_error("non-finite or unbounded Whitney position");
    AxisWeights along[3], transverse[3];
    for (int a = 0; a < 3; ++a) {
        along[a] = quadratic_axis(u[a]);
        // The face coordinate is staggered by half a cell in both transverse axes.
        transverse[a] = linear_axis(u[a] - 0.5);
    }
    MagneticStencil result;
    for (int c = 0; c < 3; ++c) {
        const int a = (c + 1) % 3, b = (c + 2) % 3;
        for (int ic = 0; ic < along[c].count; ++ic)
            for (int ia = 0; ia < transverse[a].count; ++ia)
                for (int ib = 0; ib < transverse[b].count; ++ib) {
                    int node[3];
                    node[c] = along[c].base + ic;
                    node[a] = transverse[a].base + ia;
                    node[b] = transverse[b].base + ib;
                    append_weight(result, node[0], node[1], node[2], c,
                                  along[c].weight[ic] * transverse[a].weight[ia] * transverse[b].weight[ib], domain);
                }
    }
    return result;
}

inline ElectricStencil electric_stencil(const Vector3R& position, const Domain& domain, KernelOrder order) {
    return order == KernelOrder::Linear ? electric_stencil_linear(position, domain)
                                        : electric_stencil_quadratic(position, domain);
}

inline MagneticStencil magnetic_stencil(const Vector3R& position, const Domain& domain, KernelOrder order) {
    return order == KernelOrder::Linear ? magnetic_stencil_linear(position, domain)
                                        : magnetic_stencil_quadratic(position, domain);
}

inline std::array<EdgeWeight, 12> electric_weights(const Vector3R& position, const Domain& domain) {
    const auto stencil = electric_stencil_linear(position, domain);
    std::array<EdgeWeight, 12> result;
    for (int i = 0; i < 12; ++i) result[i] = stencil.weight[i];
    return result;
}

inline Vector3R gather_electric(const Field3d& field, const Vector3R& position, const Domain& domain,
                                KernelOrder order = KernelOrder::Linear) {
    Vector3R result(0, 0, 0);
    const auto stencil = electric_stencil(position, domain, order);
    for (int i = 0; i < stencil.size; ++i) {
        const auto& w = stencil.weight[i];
        result[w.component] += w.weight * field[w.index];
    }
    return result;
}

inline Vector3R gather_magnetic(const Field3d& field, const Vector3R& position, const Domain& domain,
                                KernelOrder order = KernelOrder::Linear) {
    Vector3R result(0, 0, 0);
    const auto stencil = magnetic_stencil(position, domain, order);
    for (int i = 0; i < stencil.size; ++i) {
        const auto& w = stencil.weight[i];
        result[w.component] += w.weight * field[w.index];
    }
    return result;
}

// v_bar = alpha(B) [v^n + (q dt / 2m) E_bar].
inline Vector3R magnetic_response(const Vector3R& value, const Vector3R& b) {
    return (value + value.cross(b) + b * value.dot(b)) / (1 + b.squaredNorm());
}
}   // namespace ecsimWhitney
