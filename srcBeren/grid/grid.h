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
    // Возвращает индексы (i,j,k) узла поля в сетке (с ghost-ячейками),
    // соответствующего заданной мировой координате.
    Vector3I get_field_node_index(const Vector3R& world_coord, FieldType field, int component) const {
        const Vector3R shift = yee_shift(field, component);

        double inv_dx = 1.0 / cell_size().x();
        double inv_dy = 1.0 / cell_size().y();
        double inv_dz = 1.0 / cell_size().z();

        double x = (world_coord.x() - origin().x()) * inv_dx - shift.x();
        double y = (world_coord.y() - origin().y()) * inv_dy - shift.y();
        double z = (world_coord.z() - origin().z()) * inv_dz - shift.z();

        int i = static_cast<int>(std::floor(x)) + ghost_cells_;
        int j = static_cast<int>(std::floor(y)) + ghost_cells_;
        int k = static_cast<int>(std::floor(z)) + ghost_cells_;

        return Vector3I(i, j, k);
    }

    // Physical coordinates of a stored node, including the component's Yee shift.
    Vector3R get_node_position(int i, int j, int k, FieldType field, int component) const {
        const double x0 = origin_.x() + (i - ghost_cells_) * cell_size_.x();
        const double y0 = origin_.y() + (j - ghost_cells_) * cell_size_.y();
        const double z0 = origin_.z() + (k - ghost_cells_) * cell_size_.z();
        const Vector3R shift = yee_shift(field, component);
        return {x0 + shift.x() * cell_size_.x(), y0 + shift.y() * cell_size_.y(), z0 + shift.z() * cell_size_.z()};
    }

   private:
    Vector3R cell_size_;
    Vector3I num_cells_;   // количество ячеек без ghost
    Vector3I size_;        // num_cells + 2*ghost_cells + 1
    int ghost_cells_;
    double cell_volume_;
    Vector3R origin_;
};
