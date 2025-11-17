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
#include <sstream>
#include <string>
#include <vector>
#include "containers.h"
#include "parameters_map.h"
#include "random_generator.h"
class Bounds {
   public:
    struct BoundValues {
        BoundType x;
        BoundType y;
        BoundType z;

        BoundValues(BoundType x, BoundType y, BoundType z) : x(x), y(y), z(z) {}
    };
    // Default values is periodic boundaries
    Bounds()
        : lowerBounds(BoundType::PERIODIC, BoundType::PERIODIC,
                      BoundType::PERIODIC),
          upperBounds(BoundType::PERIODIC, BoundType::PERIODIC,
                      BoundType::PERIODIC) {}

    // Sets the lower and upper bound values
    void setBounds(const BoundValues& lower, const BoundValues& upper) {
        lowerBounds = lower;
        upperBounds = upper;
    }
    void setBounds(const ParametersMap& parameters) {
        lowerBounds.x =
            get_bound_from_str(parameters.get_string("BoundTypeX", 0));
        lowerBounds.y =
            get_bound_from_str(parameters.get_string("BoundTypeY", 0));
        lowerBounds.z =
            get_bound_from_str(parameters.get_string("BoundTypeZ", 0));
        upperBounds.x =
            get_bound_from_str(parameters.get_string("BoundTypeX", 1));
        upperBounds.y =
            get_bound_from_str(parameters.get_string("BoundTypeY", 1));
        upperBounds.z =
            get_bound_from_str(parameters.get_string("BoundTypeZ", 1));
    }

    BoundType get_bound_from_str(const std::string& bound_str) {
        if (bound_str == "PERIODIC"){
            return BoundType::PERIODIC;
        }
        else if (bound_str == "OPEN"){
            return BoundType::OPEN;
        } else if (bound_str == "OPEN_RADIUS") {
            return BoundType::OPEN_RADIUS;
        } else if (bound_str == "NEIGHBOUR") {
            return BoundType::NEIGHBOUR;
        } else {
            std::cout << "Invalid bound type" << std::endl;
            exit(1);
        }
    }

    bool check_correct_bounds(){
        if (lowerBounds.x == BoundType::PERIODIC ||
            upperBounds.x == BoundType::PERIODIC) {
            return lowerBounds.x == upperBounds.x;
        }
        if (lowerBounds.y == BoundType::PERIODIC ||
            upperBounds.y == BoundType::PERIODIC) {
            return lowerBounds.y == upperBounds.y;
        }
        if (lowerBounds.z == BoundType::PERIODIC ||
            upperBounds.z == BoundType::PERIODIC) {
            return lowerBounds.z == upperBounds.z;
        }
        if (lowerBounds.x == BoundType::OPEN_RADIUS ||
            upperBounds.x == BoundType::OPEN_RADIUS ||
            lowerBounds.y == BoundType::OPEN_RADIUS ||
            upperBounds.y == BoundType::OPEN_RADIUS) {
            bool is_correct_x = lowerBounds.x == upperBounds.x;
            bool is_correct_y = lowerBounds.y == upperBounds.y;
            return is_correct_x && is_correct_y && lowerBounds.x == upperBounds.y;
        }

            return true;
        }
    // Lower boundary conditions
    BoundValues lowerBounds;

    // Upper boundary conditions
    BoundValues upperBounds;

    bool isPeriodic(const int dim) const {
        switch (dim) {
            case X:
                return lowerBounds.x == BoundType::PERIODIC &&
                       upperBounds.x == BoundType::PERIODIC;
            case Y:
                return lowerBounds.y == BoundType::PERIODIC &&
                       upperBounds.y == BoundType::PERIODIC;
            case Z:
                return lowerBounds.z == BoundType::PERIODIC &&
                       upperBounds.z == BoundType::PERIODIC;
            default:
                std::cout << "Invalid dimensionin in check bound" << std::endl;
                return false;
        }
    }
    bool isOpenRadius() const{
        return lowerBounds.x == BoundType::OPEN_RADIUS &&
               upperBounds.x == BoundType::OPEN_RADIUS && 
               lowerBounds.y == BoundType::OPEN_RADIUS &&
               upperBounds.y == BoundType::OPEN_RADIUS;
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
    Domain(const ParametersMap& parameters, const Bounds& bound);
    Domain();
    void setDomain(const ParametersMap& parameters, const Bounds& bound);

    double3 cell_size() const { return mCellSize; }
    double cell_size(int dim) const { return mCellSize(dim); }
    double cell_volume() const { return mCellSize.x()*mCellSize.y()*mCellSize.z(); }
    int3 origin() const { return mOrigin; }
    void set_origin(const int3& newOrigin) { mOrigin = newOrigin; }
    int3 num_cells() const { return mNumCells; }
    int num_cells(const int dim) const { return mNumCells(dim); }
    int3 size() const { return mSize; }
    Bounds::BoundValues lower_bounds() const { return mBound.lowerBounds; }
    Bounds::BoundValues upper_bounds() const { return mBound.upperBounds; }
    bool is_periodic_bound(const int dim) const { return mBound.isPeriodic(dim); }
    Bounds get_bounds() const { return mBound; }

    bool is_ghost_cell(int i, int j, int k) const{
        return (i < GHOST_CELLS || i > mNumCells.x() - 1 + GHOST_CELLS ||
                j < GHOST_CELLS || j > mNumCells.y() - 1 + GHOST_CELLS ||
                k < GHOST_CELLS || k > mNumCells.z() - 1 + GHOST_CELLS);
    }

    std::tuple<bool, Axis> in_bbox_region(const double3& x) const {
        for (int i = 0; i < MAX_DIM; i++) {
            double xi = x(i) / mCellSize(i);
            if (xi <= 0 || xi >= mNumCells(i)) {
                return {false, static_cast<Axis>(i)};
            }
        }
        return {true, Axis::C};
    }

    std::tuple < bool, Axis > in_region_trap(const double3& x) const {
        for (int i = 0; i < MAX_DIM; i++) {
            double xi = x(i) / mCellSize(i);
            if (xi <= 0 || xi >= mNumCells(i)) {
                return {false, static_cast<Axis>(i)};
            }
        }
        if (true) {
            double Rx = 0.5*mNumCells.x()*mCellSize.x();
            double Ry = 0.5*mNumCells.y()*mCellSize.y();
            double R = std::min(Rx, Ry);
            double cx = x.x() - Rx;
            double cy = x.y() - Ry;
            if( cx*cx + cy*cy > R*R) {
                return {false, Axis::X};
            }
        }
        return {true, Axis::C};
    }

    // Проверка цилиндра вдоль Z (радиус по XY).
    // dim < 0  -> полная проверка; dim >= 0 -> проверка только для указанной
    // оси (радиус применяется, если dim != Z)
    inline std::tuple<bool, Axis> check_cylinder_z(const double3& x,
                                                   int dim) const {
        const bool full_check = (dim < 0);
        if (mBound.isOpenRadius() && (full_check || dim != Z)) {
            const double Rx = 0.5 * mNumCells(X) * mCellSize(X);
            const double Ry = 0.5 * mNumCells(Y) * mCellSize(Y);
            const double R = std::min(Rx, Ry);

            const double cx = x.x() - Rx;
            const double cy = x.y() - Ry;
            if (cx * cx + cy * cy >= R * R) {
                return {false, Axis::X};   // как в исходной семантике
            }
        }
        return {true, Axis::C};
    }

    // Быстрая inline-обёртка, если нужна только булева проверка
    inline bool check_cylinder_z_bool(const double3& x, int dim) const {
        return std::get<0>(check_cylinder_z(x, dim));
    }

    // Проверка по bounding-box (для всех осей или по одной)
    inline std::tuple<bool, Axis> check_bbox_dim(const double3& x,
                                                 int dim) const {
        const bool full_check = (dim < 0);

        auto check_one_dim = [&](int i) -> bool {
            const double xi = x(i) / mCellSize(i);

            if (i == X) {
                // Здесь используем условие, которое ты привёл (не PERIODIC ->
                // требуем попадания в (0, mNumCells])
                if (mBound.lowerBounds.x != BoundType::PERIODIC &&
                    (xi <= 0 || xi > mNumCells.x()))
                    return false;
            } else if (i == Y) {
                if (mBound.lowerBounds.y != BoundType::PERIODIC &&
                    (xi <= 0 || xi > mNumCells.y()))
                    return false;
            } else {   // Z
                if (mBound.lowerBounds.z != BoundType::PERIODIC &&
                    (xi <= 0 || xi > mNumCells.z()))
                    return false;
            }
            return true;
        };

        if (full_check) {
            for (int i = 0; i < MAX_DIM; ++i) {
                if (!check_one_dim(i)) {
                    return {false, static_cast<Axis>(i)};
                }
            }
            return {true, Axis::C};
        } else {
            if (!check_one_dim(dim))
                return {false, static_cast<Axis>(dim)};
            return {true, Axis::C};
        }
    }

    inline bool check_bbox_dim_bool(const double3& x, int dim) const {
        return std::get<0>(check_bbox_dim(x, dim));
    }

    // "Истинная" реализация: сначала цилиндр (если применимо), затем bbox.
    inline std::tuple<bool, Axis> in_region_impl(const double3& x,
                                                 int dim) const {
        // 1) cylinder check: если там ошибка — сразу вернём
        auto [ok_cyl, axis_cyl] = check_cylinder_z(x, dim);
        if (!ok_cyl)
            return {false, axis_cyl};

        // 2) bbox check — возвращаем его результат (он уже вернёт ось провала,
        // если есть)
        return check_bbox_dim(x, dim);
    }

    inline std::tuple<bool, Axis> in_region(const double3& x) const {
        return in_region_impl(x, -1);
    }

    inline bool in_region(const double3& x, int dim) const {
        return std::get<0>(in_region_impl(x, dim));
    }
    bool in_region_electric(int i, int j, int k, int d) const {
        bool in_region = true;
        if (mBound.lowerBounds.z == BoundType::OPEN) {
            // Ez, k=0  ==  -0.5*Dx
            if (d == Z && k == 0) {
                in_region = false;
            }
            // Ex, Ey, k=0  ==  -Dx
            if (d != Z && k <= 1)
                in_region = false;
        }
        if (mBound.upperBounds.z == BoundType::OPEN) {
            if (k >= mSize.z() - 2)
                in_region = false;
        }

        if (mBound.lowerBounds.x == BoundType::OPEN) {
            // Ex, i=0  ==  -0.5*Dx
            if (d == X && i == 0) {
                in_region = false;
            }
            // Ez, Ey, i=0  ==  -Dx
            if (d != X && i <= 1)
                in_region = false;
        }
        if (mBound.upperBounds.x == BoundType::OPEN) {
            if (i >= mSize.x() - 2)
                in_region = false;
        }
        if (mBound.lowerBounds.y == BoundType::OPEN) {
            // Ey, j=0  ==  -0.5*Dy
            if (d == Y && j == 0) {
                in_region = false;
            }
            // Ez, Ex, j=0  ==  -Dy
            if (d != Y && j <= 1)
                in_region = false;
        }
        if (mBound.upperBounds.y == BoundType::OPEN) {
            if (j >= mSize.y() - 2)
                in_region = false;
        }
        if (mBound.isOpenRadius()) {

            double Rx = 0.5*mNumCells.x()*mCellSize.x();
            double Ry = 0.5*mNumCells.y()*mCellSize.y();
            double R = std::min(Rx, Ry);
            double ix = (i-1)*mCellSize.x() - R;
            double iy = (j-1)*mCellSize.y() - R;
            if(d == X) {
                ix += 0.5*mCellSize.x();
            }
            if(d == Y) {
                iy += 0.5*mCellSize.y();
            }
            if(ix*ix + iy*iy >= R*R) {
              in_region = false;
            }
        }
        return in_region;
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
    bool in_region_magnetic(int i, int j, int k, int d) const {
        bool in_region = true;
        if (mBound.lowerBounds.z == BoundType::OPEN) {
            // Bz, k=0  ==  -Dx
            if (d == Z && k <= 1) {
                in_region = false;
            }
            // Ex, Ey, k=0  ==  -Dx
            if (k == 0) in_region = false;
        }
        if (mBound.upperBounds.z == BoundType::OPEN) {
            if (k >= mSize.z() - 2)
                in_region = false;
        }
        if (mBound.lowerBounds.x == BoundType::OPEN) {
            // Bz, k=0  ==  -Dx
            if (d == X && i <= 1) {
                in_region = false;
            }
            // Ex, Ey, k=0  ==  -Dx
            if (i == 0) in_region = false;
        }
        if (mBound.upperBounds.x == BoundType::OPEN) {
            if (i >= mSize.x() - 2)
                in_region = false;
        }
        if (mBound.lowerBounds.y == BoundType::OPEN) {
            // Bz, k=0  ==  -Dx
            if (d == Y && j <= 1) {
                in_region = false;
            }
            // Ex, Ey, k=0  ==  -Dx
            if (j == 0) in_region = false;
        }
        if (mBound.upperBounds.y == BoundType::OPEN) {
            if (j >= mSize.y() - 2)
                in_region = false;
        }

        if (mBound.isOpenRadius()) {
            double Rx = 0.5 * mNumCells.x() * mCellSize.x();
            double Ry = 0.5 * mNumCells.y() * mCellSize.y();
            double R = std::min(Rx, Ry);
            double ix = (i-1) * mCellSize.x() - R;
            double iy = (j-1) * mCellSize.y() - R;
            if (d == Y || d == Z) {
                ix += 0.5 * mCellSize.x();
            }
            if (d == X || d == Z) {
                iy += 0.5 * mCellSize.y();
            }
            if (ix * ix + iy * iy >= R * R) {
                in_region = false;
            }
        }
        return in_region;
    }

    bool in_region_magnetic(int index) const { 
        const int i = pos_vind(index, 0);
        const int j = pos_vind(index, 1);
        const int k = pos_vind(index, 2);
        const int d = pos_vind(index, 3);

        return in_region_magnetic(i, j, k, d);
    }


    bool in_region_electric(int index) const { 
        const int i = pos_vind(index, 0);
        const int j = pos_vind(index, 1);
        const int k = pos_vind(index, 2);
        const int d = pos_vind(index, 3);

        return in_region_electric(i, j, k, d);
    }

    bool in_region_density(int i, int j, int k) const {
        bool in_region = true;
        if (mBound.lowerBounds.z == BoundType::OPEN) {
            if(k == 0) in_region = false;
        }
        if (mBound.upperBounds.z == BoundType::OPEN) {
            if(k >= mNumCells.z() - 1) in_region = false;
        }
        if (mBound.lowerBounds.x == BoundType::OPEN) {
            if(i == 0) in_region = false;
        }
        if (mBound.upperBounds.x == BoundType::OPEN) {
            if(i >= mNumCells.x() - 1) in_region = false;
        }
        if (mBound.lowerBounds.y == BoundType::OPEN) {
            if(j == 0) in_region = false;
        }
        if (mBound.upperBounds.y == BoundType::OPEN) {
            if(j >= mNumCells.y() - 1) in_region = false;
        }
        if (mBound.isOpenRadius()) {
            double Rx = 0.5 * mNumCells.x() * mCellSize.x();
            double Ry = 0.5 * mNumCells.y() * mCellSize.y();
            double R = std::min(Rx, Ry);
            double ix = (i - 1) * mCellSize.x() - R;
            double iy = (j - 1) * mCellSize.y() - R;
            if (ix * ix + iy * iy >= R * R) {
                return false;
            }
        }
        return in_region;
    }

    bool in_region_density(int index) const { 
        const int i = pos_sind(index, 0);
        const int j = pos_sind(index, 1);
        const int k = pos_sind(index, 2);
        return in_region_density(i, j, k);
    }

    void make_point_periodic(double3& coord) const {
        for (int i = 0; i < MAX_DIM; i++) {
            if (mBound.isPeriodic(i)) {
                if (coord(i) < 0.) {
                    coord(i) += mNumCells(i) * mCellSize(i);
                }
                if (coord(i) >= mNumCells(i) * mCellSize(i)) {
                    coord(i) -= mNumCells(i) * mCellSize(i);
                }
            }
        }
}
    // /**
    //  * Checks if the given x coordinate is within the domain region
    //  */

    // bool in_region(const double x, const int dim) const {
    //     double xi = x / mCellSize(dim);
    //     if (xi < 0 || xi >= mNumCells(dim)) {
    //         return false;
    //     }
    //     return true;
    // }

    // /**
    //  * Checks if the given x coordinate is outside the boundary
    //  * for the given dimension dim.
    //  */
    // bool is_lost_left(const double x, const int dim) const {
    //     double xi = x / mCellSize(dim);
    //     if (xi < 0) {
    //         return true;
    //     }
    //     return false;
    // }
    // bool is_lost_right(const double x, const int dim) const {
    //     double xi = x / mCellSize(dim);
    //     if (xi >= mNumCells(dim)) {
    //         return true;
    //     }
    //     return false;
    // }

    // coordinate conversion methods
    double3 convert_global_to_local_coord(const double3& globalCoord) const {
        return globalCoord - double3(mOrigin(Dim::X) * mCellSize(Dim::X),
                                     mOrigin(Dim::Y) * mCellSize(Dim::Y),
                                     mOrigin(Dim::Z) * mCellSize(Dim::Z));
    }
    double convert_global_to_local_coord(const double globalCoord,
                                         const int dim) const {
        return globalCoord - mOrigin(dim) * mCellSize(dim);
    }
    double3 convert_local_to_global_coord(const double3& localCoord) const {
        return localCoord + double3(mOrigin(Dim::X) * mCellSize(Dim::X),
                                    mOrigin(Dim::Y) * mCellSize(Dim::Y),
                                    mOrigin(Dim::Z) * mCellSize(Dim::Z));
    }
    double convert_local_to_global_coord(const double localCoord,
                                         const int dim) const {
        return localCoord + mOrigin(dim) * mCellSize(dim);
    }
    int convert_global_to_local_index(const int indx, const int dim) const {
        return indx - mOrigin(dim);
    }

    int total_size() const { return mSize.x() * mSize.y() * mSize.z(); };
    
    // relative coordinate in cells taking into account ghost cells
    double3 get_coord_in_cell(const double3 &coord) const{
        return coord / mCellSize + double3(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);
    }
    int get_node_from_coord(const double coord, const int dim) const {
        return int(coord / mCellSize(dim) + GHOST_CELLS);
    }


    void get_interpolation_env(const double3 coord, int3& index, double3& weight, double shift) const;
    InterpolationEnvironment get_interpolation_environment(const double3 coord, double shift) const;
    double3 interpolate_fieldB(const Field3d& field, const double3& coord) ;
   private:
    double3 mCellSize;
    int3 mOrigin;
    int3 mNumCells;
    int3 mSize; // size + Ghosts
    Bounds mBound;
    };

#endif
