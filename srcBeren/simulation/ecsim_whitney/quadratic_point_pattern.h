#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "grid/domain.h"
#include "utils/sparse.h"
#include "whitney_shape.h"

namespace ecsimWhitney {

// Fixed row window for the quadratic point response on a fully periodic box.
//
// The eight base keys per cell have fixed 54-entry stencils. Their union gives
// a translation-invariant set of (row, column) offsets. We build that pattern
// once and then scatter the nine 18x18 response blocks directly into the CSR
// values instead of emitting and sorting Eigen triplets every step.
class QuadraticPointPattern {
   public:
    struct Delta {
        int dx = 0, dy = 0, dz = 0, component = 0;
    };

    // Reused between steps; one entry per cell processed in a chunk.
    struct Chunk {
        std::array<std::array<Eigen::Matrix<double, 18, 18>, 9>, 8> blocks;
        std::array<ElectricStencil, 8> stencils;
        std::array<bool, 8> used{};
    };

    bool build(const Domain& domain, Operator& matrix);
    int support() const {
        return support_;
    }
    int key_slot(int key, int p, int q) const {
        return key_slot_[static_cast<size_t>(key)][static_cast<size_t>(p * 54 + q)];
    }
    int csr_slot(int row, int offset_slot) const {
        return slot_of_offset_[static_cast<size_t>(row) * support_ + offset_slot];
    }
    Chunk& chunk(int index) {
        return chunks_[index];
    }
    const Chunk& chunk(int index) const {
        return chunks_[index];
    }
    void copy_alias_rows(Operator& matrix) const;

   private:
    static int canonical(int value, int period) {
        return ((value % period) + period) % period;
    }
    static bool delta_less(const Delta& a, const Delta& b) {
        return std::tie(a.dx, a.dy, a.dz, a.component) < std::tie(b.dx, b.dy, b.dz, b.component);
    }
    static bool delta_equal(const Delta& a, const Delta& b) {
        return a.dx == b.dx && a.dy == b.dy && a.dz == b.dz && a.component == b.component;
    }

    int rows_ = 0;
    int support_ = 0;
    std::array<std::array<int, 54 * 54>, 8> key_slot_{};
    std::vector<int> slot_of_offset_;
    std::vector<int> source_row_;
    std::vector<Chunk> chunks_;
};

inline bool QuadraticPointPattern::build(const Domain& domain, Operator& matrix) {
    for (int axis = 0; axis < 3; ++axis)
        if (!domain.is_periodic(axis))
            return false;

    const Vector3I cells = domain.num_cells();
    if (cells.x() < 3 || cells.y() < 3 || cells.z() < 3)
        return false;
    const Vector3I cell(cells.x() / 2, cells.y() / 2, cells.z() / 2);
    const Vector3I base_storage(cell.x() + GHOST_CELLS, cell.y() + GHOST_CELLS, cell.z() + GHOST_CELLS);

    std::array<std::array<Delta, 54>, 8> key_dof;
    for (int key = 0; key < 8; ++key) {
        Vector3R world;
        for (int axis = 0; axis < 3; ++axis) {
            const double fraction = (key & (1 << axis)) ? 0.75 : 0.25;
            world[axis] = (cell[axis] + fraction) * domain.cell_size(axis);
        }
        const auto stencil = electric_stencil_quadratic(world, domain);
        if (stencil.size != 54)
            return false;
        for (int p = 0; p < 54; ++p) {
            const auto [i, j, k, component] = domain.grid.pos_vind_range(stencil.weight[p].index);
            key_dof[key][p] = {i - base_storage.x(), j - base_storage.y(), k - base_storage.z(), component};
        }
    }

    std::array<std::vector<Delta>, 3> offsets;
    for (int key = 0; key < 8; ++key)
        for (int p = 0; p < 54; ++p) {
            const int row_component = key_dof[key][p].component;
            for (int q = 0; q < 54; ++q) {
                Delta delta{key_dof[key][q].dx - key_dof[key][p].dx,
                            key_dof[key][q].dy - key_dof[key][p].dy,
                            key_dof[key][q].dz - key_dof[key][p].dz,
                            key_dof[key][q].component};
                delta.dx = canonical(delta.dx, cells.x());
                delta.dy = canonical(delta.dy, cells.y());
                delta.dz = canonical(delta.dz, cells.z());
                offsets[row_component].push_back(delta);
            }
        }
    for (auto& list : offsets) {
        std::sort(list.begin(), list.end(), delta_less);
        list.erase(std::unique(list.begin(), list.end(), delta_equal), list.end());
    }
    for (int component = 0; component < 3; ++component)
        support_ = std::max(support_, static_cast<int>(offsets[component].size()));
    if (support_ == 0)
        return false;

    for (int key = 0; key < 8; ++key)
        for (int p = 0; p < 54; ++p) {
            const int row_component = key_dof[key][p].component;
            for (int q = 0; q < 54; ++q) {
                Delta delta{key_dof[key][q].dx - key_dof[key][p].dx,
                            key_dof[key][q].dy - key_dof[key][p].dy,
                            key_dof[key][q].dz - key_dof[key][p].dz,
                            key_dof[key][q].component};
                delta.dx = canonical(delta.dx, cells.x());
                delta.dy = canonical(delta.dy, cells.y());
                delta.dz = canonical(delta.dz, cells.z());
                const auto& list = offsets[row_component];
                const auto it = std::lower_bound(list.begin(), list.end(), delta, delta_less);
                if (it == list.end() || !delta_equal(*it, delta))
                    return false;
                key_slot_[key][p * 54 + q] = static_cast<int>(it - list.begin());
            }
        }

    rows_ = 3 * domain.total_size();
    slot_of_offset_.assign(static_cast<size_t>(rows_) * support_, 0);
    source_row_.assign(rows_, 0);
    std::vector<Trip> entries;
    entries.reserve(static_cast<size_t>(rows_) * support_);

    std::vector<int> columns(support_), order(support_);
    for (int row = 0; row < rows_; ++row) {
        const auto [i, j, k, component] = domain.grid.pos_vind_range(row);
        const auto& list = offsets[component];
        for (int t = 0; t < support_; ++t) {
            if (t < static_cast<int>(list.size())) {
                const auto& delta = list[t];
                const int ci = shape_index(i + delta.dx, 0, domain);
                const int cj = shape_index(j + delta.dy, 1, domain);
                const int ck = shape_index(k + delta.dz, 2, domain);
                columns[t] = domain.vind(ci, cj, ck, delta.component);
            } else {
                // Padding keeps a uniform row length; key_slot never selects it.
                columns[t] = row;
            }
        }
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            if (columns[a] != columns[b])
                return columns[a] < columns[b];
            return a < b;
        });
        for (int s = 0; s < support_; ++s) {
            slot_of_offset_[static_cast<size_t>(row) * support_ + order[s]] = s;
            entries.emplace_back(row, columns[order[s]], 0.0);
        }
        const int si = shape_index(i, 0, domain), sj = shape_index(j, 1, domain),
                  sk = shape_index(k, 2, domain);
        source_row_[row] = domain.vind(si, sj, sk, component);
    }

    matrix.resize(rows_, rows_);
    matrix.setFromTriplets(entries.begin(), entries.end());
    matrix.makeCompressed();
    if (matrix.nonZeros() != static_cast<Eigen::Index>(rows_) * support_)
        return false;
    chunks_.resize(128);
    return true;
}

inline void QuadraticPointPattern::copy_alias_rows(Operator& matrix) const {
    double* values = matrix.valuePtr();
    if (!values)
        throw std::runtime_error("quadratic Whitney pattern is not compressed");
#pragma omp parallel for schedule(static)
    for (int row = 0; row < rows_; ++row) {
        const int source = source_row_[row];
        if (source == row)
            continue;
        for (int t = 0; t < support_; ++t) {
            const int dst = slot_of_offset_[static_cast<size_t>(row) * support_ + t];
            const int src = slot_of_offset_[static_cast<size_t>(source) * support_ + t];
            values[static_cast<size_t>(row) * support_ + dst] = values[static_cast<size_t>(source) * support_ + src];
        }
    }
}

}   // namespace ecsimWhitney
