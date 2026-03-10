#ifndef WORLD_H_
#define WORLD_H_
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
#include "service.h"

// Структура для цилиндрического выреза
struct CylindricalCut {
    bool active = false;
    Vector3R center;   // центр цилиндра (обычно центр области)
    double radius = 0.0;

    bool is_inside(const Vector3R& pos) const {
        if (!active)
            return true;
        double dx = pos.x() - center.x();
        double dy = pos.y() - center.y();
        return (dx * dx + dy * dy) <= radius * radius;
    }
    // Частично или полностью вне круга
    bool is_cell_outside(double xmin, double xmax,
                                              double ymin, double ymax) const {
        // Квадрат радиуса для сравнения без извлечения корня
        double r2 = radius * radius;
        double maxDist2 = 0.0;

        // Проверяем четыре вершины
        auto checkVertex = [&](double x, double y) {
            double dx = x - center.x();
            double dy = y - center.y();
            double d2 = dx * dx + dy * dy;
            if (d2 > maxDist2)
                maxDist2 = d2;
        };

        checkVertex(xmin, ymin);
        checkVertex(xmin, ymax);
        checkVertex(xmax, ymin);
        checkVertex(xmax, ymax);

        // Если хотя бы одна вершина вне или на границе -> true
        return maxDist2 >= r2;
    }
};

enum class FieldType { ELECTRIC, MAGNETIC, DENSITY, CURRENT };
enum class Side { LOWER, UPPER };

class Bounds {
   public:
    std::array<BoundType, 3> lower;   // по индексам X, Y, Z
    std::array<BoundType, 3> upper;
    CylindricalCut cylinder;

    // Конструктор по умолчанию: все периодические
    Bounds() {
        lower.fill(BoundType::PERIODIC);
        upper.fill(BoundType::PERIODIC);
    }

    // ---- snake_case имена ----
    void set_bounds(const Bounds& bounds) {
        lower = bounds.lower;
        upper = bounds.upper;
        cylinder = bounds.cylinder;
    }

    void set_bounds(const std::array<BoundType, 3>& lower_bound,
                    const std::array<BoundType, 3>& upper_bound) {
        lower = lower_bound;   // исправлено: ранее параметр затенял поле
        upper = upper_bound;
    }

    // Загрузка из JSON
    void set_bounds(const nlohmann::json& config) {
        try {
            if (config.contains("BoundTypeX") &&
                config["BoundTypeX"].is_array() &&
                config["BoundTypeX"].size() == 2) {
                lower[0] = get_bound_from_str(
                    config["BoundTypeX"][0].get<std::string>());
                upper[0] = get_bound_from_str(
                    config["BoundTypeX"][1].get<std::string>());
            }
            if (config.contains("BoundTypeY") &&
                config["BoundTypeY"].is_array() &&
                config["BoundTypeY"].size() == 2) {
                lower[1] = get_bound_from_str(
                    config["BoundTypeY"][0].get<std::string>());
                upper[1] = get_bound_from_str(
                    config["BoundTypeY"][1].get<std::string>());
            }
            if (config.contains("BoundTypeZ") &&
                config["BoundTypeZ"].is_array() &&
                config["BoundTypeZ"].size() == 2) {
                lower[2] = get_bound_from_str(
                    config["BoundTypeZ"][0]
                        .get<std::string>());   // исправлено: индекс 2
                upper[2] = get_bound_from_str(
                    config["BoundTypeZ"][1]
                        .get<std::string>());   // исправлено: индекс 2
            }
        } catch (const nlohmann::json::exception& e) {
            std::cerr << "Error: Boundary conditions was not set" << std::endl;
            exit(-1);
        }

        // Чтение цилиндрического выреза
        if (config.contains("CylindricalCutXY") &&
            config["CylindricalCutXY"].is_object()) {
            const auto& cyl = config["CylindricalCutXY"];
            cylinder.active = true;
            if (cylinder.active) {
                cylinder.radius = cyl["radius"].get<double>();
                if (cyl.contains("center")) {
                    cylinder.center =
                        Vector3R(cyl["center"][0], cyl["center"][1], 0);
                    std::cout << "cyl: radius " << cylinder.radius
                              << ". center: " << cylinder.center.x() << ", "
                              << cylinder.center.y() << "\n";
                }
            }
        }
    }

    bool is_periodic(int dim) const {
        return lower[dim] == BoundType::PERIODIC &&
               upper[dim] == BoundType::PERIODIC;
    }

    // Проверка, является ли граница открытой (для обратной совместимости)
    bool is_open(int dim) const {
        return lower[dim] == BoundType::OPEN || upper[dim] == BoundType::OPEN;
    }

   private:
    // Статический вспомогательный метод
    static BoundType get_bound_from_str(const std::string& bound_str) {
        if (bound_str == "PERIODIC")
            return BoundType::PERIODIC;
        if (bound_str == "OPEN")
            return BoundType::OPEN;
        if (bound_str == "NEIGHBOUR")
            return BoundType::NEIGHBOUR;
        std::cerr << "Invalid bound type: " << bound_str << std::endl;
        exit(1);
    }
};

struct InterpolationEnvironment {
    int xIndex, yIndex, zIndex;
    alignas(64) double xWeight[2], yWeight[2], zWeight[2];
};

/**
 * Domain class represents the spatial domain and grid structure.
 * It contains parameters like cell size, number of cells, boundary conditions,
 * and provides utility methods for coordinate conversion between global and
 * local indices.
 */
class Domain {
   public:
    Domain(const nlohmann::json& config, const Bounds& bound);
    Domain();
    void set_domain(const nlohmann::json& config, const Bounds& bound);
    void set_domain(const nlohmann::json& config);

    Vector3R cell_size() const { return mCellSize; }
    double cell_size(int dim) const { return mCellSize[dim]; }
    double cell_volume() const { return mCellSize.x()*mCellSize.y()*mCellSize.z(); }
    Vector3I origin() const { return mOrigin; }
    void set_origin(const Vector3I& newOrigin) { mOrigin = newOrigin; }
    Vector3I num_cells() const { return mNumCells; }
    int num_cells(const int dim) const { return mNumCells[dim]; }
    Vector3I size() const { return mSize; }
    const Bounds& get_bounds() const { return bounds_; }

    bool is_periodic_bound(const int dim) const {
        return bounds_.is_periodic(dim);
    }

    bool is_ghost_cell(int i, int j, int k) const{
        return (i < GHOST_CELLS || i > mNumCells.x() - 1 + GHOST_CELLS ||
                j < GHOST_CELLS || j > mNumCells.y() - 1 + GHOST_CELLS ||
                k < GHOST_CELLS || k > mNumCells.z() - 1 + GHOST_CELLS);
    }

    bool is_cell_outside(int i, int j, int k) const {
        if (bounds_.cylinder.active)
            return bounds_.cylinder.is_cell_outside(
                       i * cell_size().x(), (i + 1) * cell_size().y(),
                       j * cell_size().y(), (j + 1) * cell_size().y()) ||
                   is_ghost_cell(i, j, k);

        return is_ghost_cell(i, j, k);
    }
    int neighbor_index_periodic(int current, int direction, int size,
                                int ghostCells) const {
        // direction: +1 или -1
        int next = current + direction;
        if (next < 0)
            next = size - 2 * ghostCells - 1 - 1;   // пример логики
        if (next >= size)
            next = 2 * ghostCells + 1;
        return next;
    }

    int neighbor_index(int current, int dim, int direction) const {
        if (bounds_.is_periodic(dim)) {
            return neighbor_index_periodic(current, direction, mSize[dim],
                                           GHOST_CELLS);
        }
        return current + direction;
    }

    // Метод для проверки цилиндра (добавлен, так как использовался в
    // in_region_impl)
    std::tuple<bool, Axis> check_cylinder(const Vector3R& x,
                                          int /*dim*/) const {
        if (bounds_.cylinder.active && !bounds_.cylinder.is_inside(x)) {
            return {false, Axis::C};   // можно указать ось, но цилиндр не
                                       // связан с конкретной осью
        }
        return {true, Axis::C};
    }

    // Проверка bounding box с учётом непериодических границ (обе стороны)
    std::tuple<bool, Axis> check_bbox_dim(const Vector3R& x, int dim) const {
        const bool full_check = (dim < 0);

        auto check_one_dim = [&](int i) -> bool {
            double xi = x[i] / mCellSize[i];
            bool lower_ok = true, upper_ok = true;

            if (bounds_.lower[i] != BoundType::PERIODIC) {
                lower_ok = (xi > 0);   // строго больше, так как на границе
                                       // может быть особый случай
            }
            if (bounds_.upper[i] != BoundType::PERIODIC) {
                upper_ok = (xi < mNumCells[i]);   // строго меньше
            }
            return lower_ok && upper_ok;
        };

        if (full_check) {
            for (int i = 0; i < MAX_DIM; ++i) {
                if (!check_one_dim(i))
                    return {false, static_cast<Axis>(i)};
            }
            return {true, Axis::C};
        } else {
            if (!check_one_dim(dim))
                return {false, static_cast<Axis>(dim)};
            return {true, Axis::C};
        }
    }

    bool check_bbox_dim_bool(const Vector3R& x, int dim) const {
        return std::get<0>(check_bbox_dim(x, dim));
    }

    // "Истинная" реализация: сначала цилиндр (если применимо), затем bbox.
    std::tuple<bool, Axis> in_region_impl(const Vector3R& x, int dim) const {
        auto [ok_cyl, axis_cyl] = check_cylinder(x, dim);
        if (!ok_cyl)
            return {false, axis_cyl};
        return check_bbox_dim(x, dim);
    }

    std::tuple<bool, Axis> in_region(const Vector3R& x) const {
        return in_region_impl(x, -1);
    }

    bool in_region(const Vector3R& x, int dim) const {
        return std::get<0>(in_region_impl(x, dim));
    }

    inline int pos_vind(int index, int n) const{
        std::vector<int> dim = {mSize.x(), mSize.y(), mSize.z(), 3};
        int capacity = 1;
        for(unsigned int i = n + 1; i < dim.size(); i++){
            capacity *= dim[i];
        }
        return (index / capacity) % dim[n];
    }
    inline int pos_sind(int index, int n) const {
        std::vector<int> dim = {mSize.x(), mSize.y(), mSize.z()};
        int capacity = 1;
        for (unsigned int i = n + 1; i < dim.size(); i++) {
            capacity *= dim[i];
        }
        return (index / capacity) % dim[n];
    }
    inline int sind(int i, int j, int k) const {
        return i * mSize.y() * mSize.z() + j * mSize.z() + k;
    };
    // index for 3D vector fields
    inline int vind(int i, int j, int k, int d, int nd = 3) const {
        return d + nd * (i * mSize.y() * mSize.z() + j * mSize.z() + k);
    };

    bool in_region_density(int index) const { 
        const int i = pos_sind(index, 0);
        const int j = pos_sind(index, 1);
        const int k = pos_sind(index, 2);
        return is_inside(i, j, k, FieldType::DENSITY, 0); //in_region_density(i, j, k);
    }

    // Проверка, находится ли узел (i,j,k) для данного поля и компоненты внутри
    // расчётной области
    bool is_inside(int i, int j, int k, FieldType field, int component) const {
        // 1. Цилиндрический вырез
        if (bounds_.cylinder.active) {
            Vector3R pos = get_idx_field_position(i, j, k, field, component);
            if (!bounds_.cylinder.is_inside(pos))
                return false;
        }

        // 2. Проверка по каждой оси
        if (!is_inside_side(i, X, Side::LOWER, field, component))
            return false;
        if (!is_inside_side(i, X, Side::UPPER, field, component))
            return false;
        if (!is_inside_side(j, Y, Side::LOWER, field, component))
            return false;
        if (!is_inside_side(j, Y, Side::UPPER, field, component))
            return false;
        if (!is_inside_side(k, Z, Side::LOWER, field, component))
            return false;
        if (!is_inside_side(k, Z, Side::UPPER, field, component))
            return false;

        return true;
    }

    // Перегрузка для линейного индекса (поле с компонентами)
    bool is_inside(int linearIndex, FieldType field) const {
        if (field == FieldType::DENSITY) {
            int i = pos_sind(linearIndex, 0);
            int j = pos_sind(linearIndex, 1);
            int k = pos_sind(linearIndex, 2);
            return is_inside(i, j, k, FieldType::DENSITY, 0);
        }
        int i = pos_vind(linearIndex, 0);
        int j = pos_vind(linearIndex, 1);
        int k = pos_vind(linearIndex, 2);
        int d = pos_vind(linearIndex, 3);
        return is_inside(i, j, k, field, d);
    }

    // Проверка одной стороны для открытой границы (OPEN)
    bool is_inside_side_open(int idx, int dim, FieldType field, int component,
                             Side side) const {
        switch (side) {
            case Side::LOWER:
                if (field == FieldType::ELECTRIC) {
                    if (component == dim) {   // нормальная компонента
                        return idx > 0;       // Ex, i=0 == -Dx/2
                    } else {                  // тангенциальные
                        return idx > 1;       // Ey, i=0 == -Dx
                    }
                } else if (field == FieldType::MAGNETIC) {
                    if (component == dim) {   // нормальная компонента
                        return idx > 1;       // Bx, i=0 == -Dx
                    } else {
                        return idx > 0;   // тангенциальные: Bx, j=0 == -Dy/2
                    }
                } else if (field == FieldType::DENSITY) {
                    return idx > 0;   // i=0 == -Dx/2
                } else {
                    return false;
                }
                break;
            case Side::UPPER:
                // i=size-1 == size + 3*Dx/2 || size + Dx
                return idx < mSize[dim] - 2;
            default:
                return true;
        }
    }
    void make_point_periodic(Vector3R& coord) const {
        for (int i = 0; i < MAX_DIM; i++) {
            if (bounds_.is_periodic(i)) {
                if (coord[i] < 0.) {
                    coord[i] += mNumCells[i] * mCellSize[i];
                }
                if (coord[i] >= mNumCells[i] * mCellSize[i]) {
                    coord[i] -= mNumCells[i] * mCellSize[i];
                }
            }
        }
}

    // coordinate conversion methods
    Vector3R convert_global_to_local_coord(const Vector3R& globalCoord) const {
        return globalCoord - Vector3R(mOrigin[Dim::X] * mCellSize[Dim::X],
                                     mOrigin[Dim::Y] * mCellSize[Dim::Y],
                                     mOrigin[Dim::Z] * mCellSize[Dim::Z]);
    }
    double convert_global_to_local_coord(const double globalCoord,
                                         const int dim) const {
        return globalCoord - mOrigin[dim] * mCellSize[dim];
    }
    Vector3R convert_local_to_global_coord(const Vector3R& localCoord) const {
        return localCoord + Vector3R(mOrigin[Dim::X] * mCellSize[Dim::X],
                                    mOrigin[Dim::Y] * mCellSize[Dim::Y],
                                    mOrigin[Dim::Z] * mCellSize[Dim::Z]);
    }
    double convert_local_to_global_coord(const double localCoord,
                                         const int dim) const {
        return localCoord + mOrigin[dim] * mCellSize[dim];
    }
    int convert_global_to_local_index(const int indx, const int dim) const {
        return indx - mOrigin[dim];
    }

    int total_size() const { return mSize.x() * mSize.y() * mSize.z(); };
    Vector3R to_cell_coordinates(const Vector3R& world_coord) const {
        return Vector3R(world_coord.x() / mCellSize.x(),
                        world_coord.y() / mCellSize.y(),
                        world_coord.z() / mCellSize.z());
    }
    Vector3I get_cell_index(const Vector3R& coord) const {
        return Vector3I{int(coord.x() / mCellSize.x() + GHOST_CELLS),
                        int(coord.y() / mCellSize.y() + GHOST_CELLS),
                        int(coord.z() / mCellSize.z() + GHOST_CELLS)};
    }

    void get_interpolation_env(const Vector3R coord, Vector3I& index, Vector3R& weight, double shift) const;
    InterpolationEnvironment get_interpolation_environment(const Vector3R coord, double shift) const;
    Vector3R interpolate_fieldB(const Field3d& field, const Vector3R& coord) ;
    bool is_inside_side(int idx, int dim, Side side, FieldType field,
                        int component) const {
        BoundType type =
            (side == Side::LOWER) ? bounds_.lower[dim] : bounds_.upper[dim];
        switch (type) {
            case BoundType::PERIODIC:
                return true;
            case BoundType::OPEN:
                return is_inside_side_open(idx, dim, field, component, side);
            case BoundType::NEIGHBOUR:
                std::cerr << "is_inside_side: not implemented for NEIGHBOUR"
                          << std::endl;
                exit(-1);
            default:
                return true;
        }
    }

   private:
    Vector3R mCellSize;
    Vector3I mOrigin;
    Vector3I mNumCells;
    Vector3I mSize; // size + Ghosts
    Bounds bounds_;
    Vector3R get_idx_field_position(int i, int j, int k, FieldType field,
                                    int component) const {
        double x = (i - GHOST_CELLS) * mCellSize.x();
        double y = (j - GHOST_CELLS) * mCellSize.y();
        double z = (k - GHOST_CELLS) * mCellSize.z();

        if (field == FieldType::ELECTRIC) {
            if (component == X)
                x += 0.5 * mCellSize.x();
            else if (component == Y)
                y += 0.5 * mCellSize.y();
            else if (component == Z)
                z += 0.5 * mCellSize.z();
        } else if (field == FieldType::MAGNETIC) {
            if (component == Y || component == Z)
                x += 0.5 * mCellSize.x();
            if (component == X || component == Z)
                y += 0.5 * mCellSize.y();
            if (component == X || component == Y)
                z += 0.5 * mCellSize.z();
        }
        // for DENSITY we do not need to shift the position 
        return {x, y, z};
    }
};

#endif
