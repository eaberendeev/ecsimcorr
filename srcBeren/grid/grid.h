#pragma once

#include <array>
#include <cmath>

#include "core/indexing.h"
#include "grid/yee.h"

class Grid {
   public:
    Grid()
        : cell_size_(0, 0, 0),
          num_cells_(0, 0, 0),
          size_(0, 0, 0),
          ghost_cells_(0),
          cell_volume_(0.0),
          origin_(0, 0, 0) {
    }

    Grid(const Vector3R& cell_size, const Vector3I& num_cells, int ghost_cells) {
        init(cell_size, num_cells, ghost_cells);
    }

    void init(const Vector3R& cell_size, const Vector3I& num_cells, int ghost_cells) {
        cell_size_ = cell_size;
        num_cells_ = num_cells;
        ghost_cells_ = ghost_cells;
        size_ = Vector3I(num_cells.x() + 2 * ghost_cells + 1, num_cells.y() + 2 * ghost_cells + 1,
                         num_cells.z() + 2 * ghost_cells + 1);
        cell_volume_ = cell_size.x() * cell_size.y() * cell_size.z();
        origin_ = Vector3R(0, 0, 0);
    }
    const Vector3R& cell_size() const {
        return cell_size_;
    }
    double cell_size(int dim) const {
        return cell_size_[dim];
    }
    const Vector3I& num_cells() const {
        return num_cells_;
    }
    int num_cells(int dim) const {
        return num_cells_[dim];
    }
    const Vector3I& size() const {
        return size_;
    }
    int size(int dim) const {
        return size_[dim];
    }
    int ghost_cells() const {
        return ghost_cells_;
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
    inline std::array<int, 4> pos_vind_range(int index) const {
        return {pos_vind(index, 0), pos_vind(index, 1), pos_vind(index, 2), pos_vind(index, 3)};
    }
    inline int pos_sind(int index, int axis) const {
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
    Vector3R to_cell_coordinates(const Vector3R& world_coord) const {
        return Vector3R(world_coord.x() / cell_size_.x(), world_coord.y() / cell_size_.y(),
                        world_coord.z() / cell_size_.z());
    }
    Vector3I get_cell_index(const Vector3R& coord) const {
        return Vector3I{int(coord.x() / cell_size_.x() + ghost_cells_), int(coord.y() / cell_size_.y() + ghost_cells_),
                        int(coord.z() / cell_size_.z() + ghost_cells_)};
    }
    // Logical node zero is the unshifted node at origin; storage starts at -ghost_cells.
    Vector3I to_storage_index(const Vector3I& node) const {
        return node + Vector3I(ghost_cells_);
    }
    Vector3I to_node_index(const Vector3I& storage) const {
        return storage - Vector3I(ghost_cells_);
    }

    // Lower logical node on this component's Yee lattice. This does not select
    // the support of a particle shape, which may require a different stencil.
    Vector3I field_node_index(const Vector3R& world_coord, FieldType field, int component) const {
        const Vector3R shift = yee_shift(field, component);
        const double x = (world_coord.x() - origin().x()) * (1.0 / cell_size().x()) - shift.x();
        const double y = (world_coord.y() - origin().y()) * (1.0 / cell_size().y()) - shift.y();
        const double z = (world_coord.z() - origin().z()) * (1.0 / cell_size().z()) - shift.z();
        return {static_cast<int>(std::floor(x)), static_cast<int>(std::floor(y)), static_cast<int>(std::floor(z))};
    }

    Vector3R node_position(const Vector3I& node, FieldType field, int component) const {
        const double x0 = origin_.x() + node.x() * cell_size_.x();
        const double y0 = origin_.y() + node.y() * cell_size_.y();
        const double z0 = origin_.z() + node.z() * cell_size_.z();
        const Vector3R shift = yee_shift(field, component);
        return {x0 + shift.x() * cell_size_.x(), y0 + shift.y() * cell_size_.y(), z0 + shift.z() * cell_size_.z()};
    }

    // Compatibility API: these indices still refer to storage, including ghosts.
    Vector3I get_field_node_index(const Vector3R& world_coord, FieldType field, int component) const {
        return to_storage_index(field_node_index(world_coord, field, component));
    }
    Vector3R get_node_position(int i, int j, int k, FieldType field, int component) const {
        return node_position(to_node_index({i, j, k}), field, component);
    }

   private:
    Vector3R cell_size_;
    Vector3I num_cells_;   // количество ячеек без ghost
    Vector3I size_;        // num_cells + 2*ghost_cells + 1
    int ghost_cells_;
    double cell_volume_;
    Vector3R origin_;
};
