#pragma once

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "containers.h"
#include "random_generator.h"
#include "util.h"

enum class FieldType { ELECTRIC, MAGNETIC, DENSITY, CURRENT };

enum class Face { XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, CYLINDER };
// enum class FieldType { ELECTRIC, MAGNETIC, DENSITY, CURRENT };

const std::vector<Face> ALL_FACES = {Face::XMIN,    Face::XMAX, Face::YMIN,
                                     Face::YMAX,    Face::ZMIN, Face::ZMAX,
                                     Face::CYLINDER};

// Смещения для узлов полей относительно центра ячейки
inline Vector3R electric_shift(int component) {
    // Ex сдвинута на +dx/2 по X, Ey на +dy/2 по Y, Ez на +dz/2 по Z
    Vector3R s{0, 0, 0};
    if (component == 0)
        s.x() = 0.5;
    else if (component == 1)
        s.y() = 0.5;
    else if (component == 2)
        s.z() = 0.5;
    return s;
}

inline Vector3R magnetic_shift(int component) {
    // Bx сдвинута на +dy/2 и +dz/2, By на +dx/2 и +dz/2, Bz на +dx/2 и +dy/2
    Vector3R s{0, 0, 0};
    if (component == 0) {
        s.y() = 0.5;
        s.z() = 0.5;
    } else if (component == 1) {
        s.x() = 0.5;
        s.z() = 0.5;
    } else if (component == 2) {
        s.x() = 0.5;
        s.y() = 0.5;
    }
    return s;
}

inline Vector3R density_shift() {
    return {0, 0, 0};   // плотность в центре ячейки
}

class Grid {
   public:
    Grid()
        : cell_size_(0, 0, 0),
          num_cells_(0, 0, 0),
          size_(0, 0, 0),
          ghost_cells_(0),
          cell_volume_(0.0),
          origin_(0, 0, 0) {}

    Grid(const Vector3R& cell_size, const Vector3I& num_cells,
         int ghost_cells) {
        init(cell_size, num_cells, ghost_cells);
    }

    void init(const Vector3R& cell_size, const Vector3I& num_cells,
              int ghost_cells) {
        cell_size_ = cell_size;
        num_cells_ = num_cells;
        ghost_cells_ = ghost_cells;
        size_ = Vector3I(num_cells.x() + 2 * ghost_cells + 1,
                         num_cells.y() + 2 * ghost_cells + 1,
                         num_cells.z() + 2 * ghost_cells + 1);
        cell_volume_ = cell_size.x() * cell_size.y() * cell_size.z();
        dims_ = {size_.x(), size_.y(), size_.z(), 3};
        origin_ = Vector3R(0, 0, 0);
    }
    const Vector3R& cell_size() const { return cell_size_; }
    double cell_size(int dim) const { return cell_size_[dim]; }
    const Vector3I& num_cells() const { return num_cells_; }
    int num_cells(int dim) const { return num_cells_[dim]; }
    const Vector3I& size() const { return size_; }
    int size(int dim) const { return size_[dim]; }
    int ghost_cells() const { return ghost_cells_; }
    double cell_volume() const { return cell_volume_; }
    const Vector3R& origin() const { return origin_; }

    inline int pos_vind(int index, int n) const {
        int capacity = 1;
        for (unsigned int i = n + 1; i < dims_.size(); i++) {
            capacity *= dims_[i];
        }
        return (index / capacity) % dims_[n];
    }
    inline int pos_sind(int index, int n) const {
        int capacity = 1;
        for (unsigned int i = n + 1; i < dims_.size() - 1; i++) {
            capacity *= dims_[i];
        }
        return (index / capacity) % dims_[n];
    }
    inline int sind(int i, int j, int k) const {
        return i * size_.y() * size_.z() + j * size_.z() + k;
    };
    // index for 3D vector fields
    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return d + nd * (i * size_.y() * size_.z() + j * size_.z() + k);
    };

    int total_size() const { return size_.x() * size_.y() * size_.z(); };
    Vector3R to_cell_coordinates(const Vector3R& world_coord) const {
        return Vector3R(world_coord.x() / cell_size_.x(),
                        world_coord.y() / cell_size_.y(),
                        world_coord.z() / cell_size_.z());
    }
    Vector3I get_cell_index(const Vector3R& coord) const {
        return Vector3I{int(coord.x() / cell_size_.x() + ghost_cells_),
                        int(coord.y() / cell_size_.y() + ghost_cells_),
                        int(coord.z() / cell_size_.z() + ghost_cells_)};
    }

   private:
    Vector3R cell_size_;
    Vector3I num_cells_;   // количество ячеек без ghost
    Vector3I size_;        // num_cells + 2*ghost_cells + 1
    std::vector<int> dims_;
    int ghost_cells_;
    double cell_volume_;
    Vector3R origin_;
};

struct Cylinder {
    Vector3R center;
    double radius;
};
// -----------------------------------------------------------------------------
// Геометрия домена: прямоугольный параллелепипед + цилиндр (внутренняя область)
// -----------------------------------------------------------------------------
struct Geometry {
    Vector3R box_min, box_max;   // границы прямоугольной области
    bool use_cylinder = false;   // активен ли цилиндр
    Cylinder cylinder;
    Vector3R cyl_center;
    double cyl_radius;

    Geometry()
        : box_min(0, 0, 0),
          box_max(0, 0, 0),
          use_cylinder(false),
          cyl_center(0, 0, 0),
          cyl_radius(0.0) {}

    Geometry(const Vector3R& bmin, const Vector3R& bmax,
             const bool cyl = false) {
        init(bmin, bmax, cyl);
    }
    void init(const Vector3R& bmin, const Vector3R& bmax,
              const bool cyl = false) {
        box_min = bmin;
        box_max = bmax;
        use_cylinder = cyl;
    }
    bool in_cylinder(const Vector3R& p, double eps = 0.0) const {
        double dx = p.x() - cyl_center.x();
        double dy = p.y() - cyl_center.y();
        return (dx * dx + dy * dy <= cyl_radius * cyl_radius - eps);
    }
    // Проверка, находится ли точка внутри области (box ∩ cylinder)
    bool contains(const Vector3R& p, double eps = 0.0) const {
        // Проверка прямоугольного параллелепипеда
        if (p.x() < box_min.x() + eps || p.x() >= box_max.x() - eps)
            return false;
        if (p.y() < box_min.y() + eps || p.y() >= box_max.y() - eps)
            return false;
        if (p.z() < box_min.z() + eps || p.z() >= box_max.z() - eps)
            return false;

        // Если цилиндр не используется, точка внутри
        if (!use_cylinder)
            return true;

        // Проверка цилиндра
        if (!in_cylinder(p, eps))
            return false;

        return true;
    }

    // Функция, возвращающая true, если точка p находится вне указанной грани
    bool is_outside_face(const Face face, const Vector3R& p,
                         double eps = 0.0) const {
        switch (face) {
            case Face::XMIN:
                return p.x() < box_min.x() + eps;
            case Face::XMAX:
                return p.x() >= box_max.x() - eps;
            case Face::YMIN:
                return p.y() < box_min.y() + eps;
            case Face::YMAX:
                return p.y() >= box_max.y() - eps;
            case Face::ZMIN:
                return p.z() < box_min.z() + eps;
            case Face::ZMAX:
                return p.z() >= box_max.z() - eps;
            case Face::CYLINDER:
                if (!use_cylinder) {
                    return false;
                }
                return !in_cylinder(p, eps);
            default:
                return false;   // на случай добавления новых граней
        }
    }
};

class Domain {
   public:
    Geometry geom;
    Grid grid;

    Vector3R cell_size() const { return grid.cell_size(); }
    double cell_size(int dim) const { return grid.cell_size(dim); }
    double cell_volume() const { return grid.cell_volume(); }
    Vector3I num_cells() const { return grid.num_cells(); }
    int num_cells(const int dim) const { return grid.num_cells(dim); }
    Vector3I size() const { return grid.size(); }
    int total_size() const {
        return grid.size().x() * grid.size().y() * grid.size().z();
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
    void init(const Vector3I& num_cells, const Vector3R& cell_size,
              int ghost = 1) {
        geom.init(Vector3R(0, 0, 0), Vector3R(num_cells.x() * cell_size.x(),
                                              num_cells.y() * cell_size.y(),
                                              num_cells.z() * cell_size.z()));
        grid.init(cell_size, num_cells, ghost);
    }

    void init_from_json(const nlohmann::json& config) {
        auto cell_size = Vector3R(get_checked<double>(config, "Dx"),
                                  get_checked<double>(config, "Dy"),
                                  get_checked<double>(config, "Dz"));

        auto num_cells = Vector3I(get_checked<int>(config, "NumCellsX"),
                                  get_checked<int>(config, "NumCellsY"),
                                  get_checked<int>(config, "NumCellsZ"));

        init(num_cells, cell_size);
        if (config.contains("CylinderDomain") &&
            config["CylinderDomain"].is_object()) {
            const auto& cyl = config["CylinderDomain"];
            auto radius = cyl["radius"].get<double>();
            if (cyl.contains("center")) {
                auto center = Vector3R(cyl["center"][0], cyl["center"][1], 0);
                std::cout << "cyl: radius " << radius
                          << ". center: " << center.x() << ", " << center.y()
                          << "\n";
                set_cylinder(center, radius);
            }
        }
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
        int i = grid.pos_vind(linear_index, 0);
        int j = grid.pos_vind(linear_index, 1);
        int k = grid.pos_vind(linear_index, 2);
        int component = field != FieldType::DENSITY ? grid.pos_vind(linear_index, 3) : 0;
        return is_inside_node(i, j, k, field, component);
    }
    // Проверка принадлежности узла поля (с учётом сдвигов Yee)
    bool is_inside_node(int i, int j, int k, FieldType field,
                        int component) const {
        // Получаем координату узла
        Vector3R pos = get_node_position(i, j, k, field, component);
        const double eps = 1.e-12;
        return contains(pos, eps);
    }

    // Получить позицию узла поля по индексам сетки (с учётом ghost)
    Vector3R get_node_position(int i, int j, int k, FieldType field,
                               int component) const {
        // Координаты левого нижнего угла ячейки (без учёта сдвига)
        double x0 =
            grid.origin().x() + (i - grid.ghost_cells()) * grid.cell_size().x();
        double y0 =
            grid.origin().y() + (j - grid.ghost_cells()) * grid.cell_size().y();
        double z0 =
            grid.origin().z() + (k - grid.ghost_cells()) * grid.cell_size().z();

        Vector3R shift;
        if (field == FieldType::ELECTRIC)
            shift = electric_shift(component);
        else if (field == FieldType::MAGNETIC)
            shift = magnetic_shift(component);
        else
            shift = density_shift();

        return {x0 + shift.x() * grid.cell_size().x(),
                y0 + shift.y() * grid.cell_size().y(),
                z0 + shift.z() * grid.cell_size().z()};
    }
};
