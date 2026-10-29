#pragma once

#include "grid/Geometry.h"
#include "grid/Grid.h"

class Domain {
   public:
    Geometry geom;
    Grid grid;

    Vector3R cell_size() const {
        return grid.cell_size();
    }
    double cell_size(int dim) const {
        return grid.cell_size(dim);
    }
    double cell_volume() const {
        return grid.cell_volume();
    }
    Vector3I num_cells() const {
        return grid.num_cells();
    }
    int num_cells(const int dim) const {
        return grid.num_cells(dim);
    }
    Vector3I size() const {
        return grid.size();
    }
    int total_size() const {
        return grid.total_size();
    };
    inline int pos_vind(int index, int n) const {
        return grid.pos_vind(index, n);
    }

    inline int pos_sind(int index, int n) const {
        return grid.pos_sind(index, n);
    }
    inline int sind(int i, int j, int k) const {
        return grid.sind(i, j, k);
    };
    // index for 3D vector fields
    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return grid.vind(i, j, k, d, nd);
    };

    Vector3I get_cell_index(const Vector3R& coord) const {
        return grid.get_cell_index(coord);
    };
    Vector3R to_cell_coordinates(const Vector3R& world_coord) const {
        return grid.to_cell_coordinates(world_coord);
    }
    // Конструктор
    Domain() = default;

    // Инициализация прямоугольной области и сетки
    void init(const Vector3I& num_cells, const Vector3R& cell_size, int ghost = 1) {
        geom.init(Vector3R(0, 0, 0), Vector3R(num_cells.x() * cell_size.x(), num_cells.y() * cell_size.y(),
                                              num_cells.z() * cell_size.z()));
        grid.init(cell_size, num_cells, ghost);
    }

    // Активация цилиндра
    void set_cylinder(const Vector3R& center, double radius) {
        geom.use_cylinder = true;
        geom.cyl_center = center;
        geom.cyl_radius = radius;
    }

    // Проверка принадлежности точки
    bool contains(const Vector3R& p, double eps = 0.0) const {
        return geom.contains(p, eps);
    }

    bool is_inside_node(int linear_index, FieldType field) const {
        const bool scalar = field == FieldType::DENSITY;
        const int i = scalar ? grid.pos_sind(linear_index, 0) : grid.pos_vind(linear_index, 0);
        const int j = scalar ? grid.pos_sind(linear_index, 1) : grid.pos_vind(linear_index, 1);
        const int k = scalar ? grid.pos_sind(linear_index, 2) : grid.pos_vind(linear_index, 2);
        const int component = scalar ? 0 : grid.pos_vind(linear_index, 3);
        return is_inside_node(i, j, k, field, component);
    }
    // Проверка принадлежности узла поля (с учётом сдвигов Yee)
    bool is_inside_node(int i, int j, int k, FieldType field, int component) const {
        // Получаем координату узла
        Vector3R pos = get_node_position(i, j, k, field, component);
        const double eps = 1.e-12;
        return contains(pos, eps);
    }
    bool is_inside_node_periodic(int i, int j, int k, FieldType field, int component, bool periodic[3]) const {
        Vector3R pos = get_node_position(i, j, k, field, component);
        return geom.contains(pos, periodic, 1e-12);
    }

    // Получить позицию узла поля по индексам сетки (с учётом ghost)
    Vector3R get_node_position(int i, int j, int k, FieldType field, int component) const {
        return grid.get_node_position(i, j, k, field, component);
    }
};
