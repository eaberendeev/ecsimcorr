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

    // Lower boundary conditions
    BoundValues lowerBounds;

    // Upper boundary conditions
    BoundValues upperBounds;
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

    /**
     * Checks if the given 3D point x is within the domain region.
     * Converts x to cell indices by dividing by cell size,
     * then checks if each index is within the number of cells in each
     * dimension.
     */
    bool in_region(const double3& x) const {
        for (int i = 0; i < MAX_DIM; i++) {
            double xi = x(i) / mCellSize(i);
            if (xi < 0 || xi >= mNumCells(i)) {
                return false;
            }
        }
        return true;
    }
    /**
     * Checks if the given x coordinate is within the domain region
     */

    bool in_region(const double x, const int dim) const {
        double xi = x / mCellSize(dim);
        if (xi < 0 || xi >= mNumCells(dim)) {
            return false;
        }
        return true;
    }

    /**
     * Checks if the given x coordinate is outside the boundary
     * for the given dimension dim.
     */
    bool is_lost_left(const double x, const int dim) const {
        double xi = x / mCellSize(dim);
        if (xi < 0) {
            return true;
        }
        return false;
    }
    bool is_lost_right(const double x, const int dim) const {
        double xi = x / mCellSize(dim);
        if (xi >= mNumCells(dim)) {
            return true;
        }
        return false;
    }

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
    int3 mSize;
    Bounds mBound;
};

#endif
