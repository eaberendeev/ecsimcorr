#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "core/indexing.h"
#include "grid/yee.h"

class Grid {
   public:
    Grid() : cell_size_(0, 0, 0), num_cells_(0, 0, 0), size_(0, 0, 0), cell_volume_(0.0), origin_(0, 0, 0) {
    }

    Grid(const Vector3R& cell_size, const Vector3I& num_cells) {
        init(cell_size, num_cells);
    }

    void init(const Vector3R& cell_size, const Vector3I& num_cells) {
        Vector3I storage_size;
        int64_t node_count = 1;
        for (Axis axis : {Axis::X, Axis::Y, Axis::Z}) {
            if (!std::isfinite(cell_size[axis]) || cell_size[axis] <= 0 || num_cells[axis] <= 0 ||
                !std::isfinite(num_cells[axis] * cell_size[axis]))
                throw std::invalid_argument("Grid: cell sizes and cell counts must be finite and positive");
            const int64_t extent = int64_t(num_cells[axis]) + 2LL * GHOST_CELLS + 1;
            // Grid strides and scalar/vector storage capacities use int indices.
            if (extent > std::numeric_limits<int>::max() || node_count > std::numeric_limits<int>::max() / 3 / extent)
                throw std::invalid_argument("Grid: scalar/vector storage exceeds the int index range");
            storage_size[axis] = static_cast<int>(extent);
            node_count *= extent;
        }
        const double volume = cell_size.x() * cell_size.y() * cell_size.z();
        if (!std::isfinite(volume) || volume <= 0)
            throw std::invalid_argument("Grid: cell volume must be finite and positive");
        cell_size_ = cell_size;
        num_cells_ = num_cells;
        size_ = storage_size;
        cell_volume_ = volume;
        origin_ = Vector3R(0, 0, 0);
    }

    const Vector3R& cell_size() const {
        return cell_size_;
    }

    double cell_size(int dim) const {
        return cell_size_[dim];
    }

    double cell_size(Axis axis) const {
        return cell_size_[axis];
    }

    const Vector3I& num_cells() const {
        return num_cells_;
    }

    int num_cells(int dim) const {
        return num_cells_[dim];
    }

    int num_cells(Axis axis) const {
        return num_cells_[axis];
    }

    const Vector3I& size() const {
        return size_;
    }

    int size(int dim) const {
        return size_[dim];
    }

    int size(Axis axis) const {
        return size_[axis];
    }

    double cell_volume() const {
        return cell_volume_;
    }

    const Vector3R& origin() const {
        return origin_;
    }

    inline int pos_vind(int index, int axis) const {
        return index_coordinate<4>(index, axis, std::array<int64_t, 4>{size_.x(), size_.y(), size_.z(), 3});
    }

    inline int pos_vind(int index, Axis axis) const {
        return index_coordinate<4>(index, axis, std::array<int64_t, 4>{size_.x(), size_.y(), size_.z(), 3});
    }

    inline std::array<int, 4> pos_vind_range(int index) const {
        return {pos_vind(index, Axis::X), pos_vind(index, Axis::Y), pos_vind(index, Axis::Z),
                pos_vind(index, Axis::Comp)};
    }

    inline int pos_sind(int index, int axis) const {
        return index_coordinate<3>(index, axis, std::array<int64_t, 3>{size_.x(), size_.y(), size_.z()});
    }

    inline int pos_sind(int index, Axis axis) const {
        return index_coordinate<3>(index, axis, std::array<int64_t, 3>{size_.x(), size_.y(), size_.z()});
    }

    inline int sind(int i, int j, int k) const {
        return ind(i, j, k, 0, size_.x(), size_.y(), size_.z(), 1);
    }

    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return ind(i, j, k, d, size_.x(), size_.y(), size_.z(), nd);
    }

    int total_size() const {
        return size_.x() * size_.y() * size_.z();
    }

    // Unique physical owner of a periodic storage index (including ghost nodes).
    int periodic_node_owner(int index, int axis) const {
        const int node = index - GHOST_CELLS;
        const int n = num_cells_[axis];
        if (node >= 0 && node < n)
            return index;
        const int wrapped = node % n;
        return (wrapped < 0 ? wrapped + n : wrapped) + GHOST_CELLS;
    }

    // Unshifted logical coordinates in cell units; ghost storage is not included.
    Vector3R to_cell_coordinates(const Vector3R& world_coord) const {
        return Vector3R(world_coord.x() / cell_size_.x(), world_coord.y() / cell_size_.y(),
                        world_coord.z() / cell_size_.z());
    }

    // Storage bin containing a world position, including halo bins.
    Vector3I storage_cell_index(const Vector3R& coord) const {
        // Floor before applying the storage offset, including below the lower halo.
        return to_storage_index(checked_floor_index(to_cell_coordinates(coord)));
    }

    bool contains_cell_index(const Vector3I& storage) const {
        for (int axis = 0; axis < 3; ++axis)
            if (storage[axis] < 0 || storage[axis] >= size_[axis])
                return false;
        return true;
    }

    // Lower logical node on this component's Yee lattice. This does not select
    // the support of a particle shape, which may require a different stencil.
    Vector3I field_node_index(const Vector3R& world_coord, FieldType field, int component) const {
        const Vector3R shift = yee_shift(field, component);
        const Vector3R normalized = to_cell_coordinates(world_coord) - shift;
        return checked_floor_index(normalized);
    }

    Vector3R node_position(const Vector3I& node, FieldType field, int component) const {
        const double x0 = origin_.x() + node.x() * cell_size_.x();
        const double y0 = origin_.y() + node.y() * cell_size_.y();
        const double z0 = origin_.z() + node.z() * cell_size_.z();
        const Vector3R shift = yee_shift(field, component);
        return {x0 + shift.x() * cell_size_.x(), y0 + shift.y() * cell_size_.y(), z0 + shift.z() * cell_size_.z()};
    }

    // Explicit storage-index API, including ghost offsets.
    Vector3I storage_field_node_index(const Vector3R& world_coord, FieldType field, int component) const {
        return to_storage_index(field_node_index(world_coord, field, component));
    }

    Vector3R storage_node_position(int i, int j, int k, FieldType field, int component) const {
        return node_position(to_node_index({i, j, k}), field, component);
    }

   private:
    Vector3R cell_size_;
    Vector3I num_cells_;   // количество ячеек без ghost
    Vector3I size_;        // num_cells + 2*GHOST_CELLS + 1
    double cell_volume_;
    Vector3R origin_;
};
