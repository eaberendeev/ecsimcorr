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

#include "const.h"
#include "containers.h"
#include "defines.h"
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


    void get_interpolation_env(const double3 coord, int3& index, double3& weight, double shift = 0.0) const {
        
        double3 coordInCell = coord / mCellSize + double3(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);

        coordInCell -= double3(shift,shift,shift);

        index.x() = static_cast<int>(coordInCell.x());
        index.y() = static_cast<int>(coordInCell.y());
        index.z() = static_cast<int>(coordInCell.z());

        weight.x() = 1 - (index.x() - coordInCell.x());
        weight.y() = 1 - (index.y() - coordInCell.y());
        weight.z() = 1 - (index.z() - coordInCell.z());
    }
    InterpolationEnvironment get_interpolation_environment(const double3 coord, double shift = 0.0) const {
        double3 coordInCell =
            coord / mCellSize + double3(GHOST_CELLS, GHOST_CELLS, GHOST_CELLS);

        coordInCell -= double3(shift, shift, shift);
        
        InterpolationEnvironment env;
        env.xIndex = static_cast<int>(coordInCell.x());
        env.yIndex = static_cast<int>(coordInCell.y());
        env.zIndex = static_cast<int>(coordInCell.z());

        env.xWeight[0] = 1 - (env.xIndex - coordInCell.x());
        env.yWeight[0] = 1 - (env.yIndex - coordInCell.y());
        env.zWeight[0] = 1 - (env.zIndex - coordInCell.z());

        env.xWeight[1] = 1 - env.xWeight[0];
        env.yWeight[1] = 1 - env.yWeight[0];
        env.zWeight[1] = 1 - env.zWeight[0];
        return env;
    }
    double3 interpolate_fieldB(const Field3d& field, const double3& coord) ;
   private:
    double3 mCellSize;
    int3 mOrigin;
    int3 mNumCells;
    int3 mSize;
    Bounds mBound;
};

struct Region {
    double3 cellSize;
    int3 cellsShift;
    double origin;
    int3 numCells;
    int3 numNodes;
    int3 dampCells[2];
    int3 boundType[2];
    bool in_region(double x) const {
        if (x < origin || x >= origin + numCells.x() * cellSize.x())
            return false;

        return true;
    }
    double3 get_coord_loc(const double3& POS) const {
        double3 POS_loc = POS;
        POS_loc.x() -= origin;
        return POS_loc;
    }
    double3 get_coord_glob(const double3& POS) const {
        double3 POS_glob = POS;
        POS_glob.x() += origin;
        return POS_glob;
    }
    //int get_index_loc(int indx) const { return indx - round(origin / Dx); }
    Region(const ParametersMap& parameters);
    int total_size() const {
        return numNodes.x() * numNodes.y() * numNodes.z();
    };
};

Region split_region(const Region& regionGlob, int rank, int splitSize);

struct World {
    World(const Region& regionGlob, const Region& regionSplit)
        : regionGlob(regionGlob), region(regionSplit){};
    Region regionGlob;
    Region region;
    ~World() {}
};

#endif
