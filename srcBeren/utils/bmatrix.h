#pragma once

#include <assert.h>

#include <unordered_map>
#include <vector>

#include "containers.h"
#include "indexing.h"
#include "sparse.h"
#include "util.h"

typedef std::unordered_map<int, double> IndexMap;

namespace BlockDims {
// Dimension constants for different components
namespace X {
constexpr int SIZE_I = 3, SIZE_J = 2, SIZE_K = 2;
}
namespace Y {
constexpr int SIZE_I = 2, SIZE_J = 3, SIZE_K = 2;
}
namespace Z {
constexpr int SIZE_I = 2, SIZE_J = 2, SIZE_K = 3;
}

// Improved indexing functions with named constants
inline constexpr int indX(int i, int j, int k) {
    return i * (X::SIZE_J * X::SIZE_K) + j * X::SIZE_K + k;
}

inline constexpr int indY(int i, int j, int k) {
    return i * (Y::SIZE_J * Y::SIZE_K) + j * Y::SIZE_K + k;
}

inline constexpr int indZ(int i, int j, int k) {
    return i * (Z::SIZE_J * Z::SIZE_K) + j * Z::SIZE_K + k;
}

// Common constants
constexpr int BLOCK_SIZE = 12;
constexpr int DIRECTIONS = 9;
constexpr int BLOCK_CAPACITY = BLOCK_SIZE * BLOCK_SIZE * DIRECTIONS;
}   // namespace BlockDims

namespace BlockDimsNGP {
// Dimension constants for different components
namespace X {
constexpr int SIZE_I = 1, SIZE_J = 2, SIZE_K = 2;
}
namespace Y {
constexpr int SIZE_I = 2, SIZE_J = 1, SIZE_K = 2;
}
namespace Z {
constexpr int SIZE_I = 2, SIZE_J = 2, SIZE_K = 1;
}

// Improved indexing functions with named constants
inline constexpr int indX(int i, int j, int k) {
    return i * (X::SIZE_J * X::SIZE_K) + j * X::SIZE_K + k;
}

inline constexpr int indY(int i, int j, int k) {
    return i * (Y::SIZE_J * Y::SIZE_K) + j * Y::SIZE_K + k;
}

inline constexpr int indZ(int i, int j, int k) {
    return i * (Z::SIZE_J * Z::SIZE_K) + j * Z::SIZE_K + k;
}

// Common constants
constexpr int BLOCK_SIZE = 4;
constexpr int DIRECTIONS = 9;
constexpr int BLOCK_CAPACITY = BLOCK_SIZE * BLOCK_SIZE * DIRECTIONS;
}   // namespace BlockDimsNGP

struct XIndexer : Indexer<3, 2, 2> {   // X: 3 точки по x, 2 по y, 2 по z
    static constexpr int dir = 0;
    static constexpr int offset_x = -1;   // Смещение по x
    static constexpr int offset_y = 0;
    static constexpr int offset_z = 0;
};
struct YIndexer : Indexer<2, 3, 2> {
    static constexpr int dir = 1;
    static constexpr int offset_x = 0;
    static constexpr int offset_y = -1;
    static constexpr int offset_z = 0;
};
struct ZIndexer : Indexer<2, 2, 3> {
    static constexpr int dir = 2;
    static constexpr int offset_x = 0;
    static constexpr int offset_y = 0;
    static constexpr int offset_z = -1;
};

struct XIndexerNGP : Indexer<1, 2, 2> {   // X: 1 точки по x, 2 по y, 2 по z
    static constexpr int dir = 0;
    static constexpr int offset_x = 0;   // Смещение по x
    static constexpr int offset_y = 0;
    static constexpr int offset_z = 0;
};
struct YIndexerNGP : Indexer<2, 1, 2> {
    static constexpr int dir = 1;
    static constexpr int offset_x = 0;
    static constexpr int offset_y = 0;
    static constexpr int offset_z = 0;
};
struct ZIndexerNGP : Indexer<2, 2, 1> {
    static constexpr int dir = 2;
    static constexpr int offset_x = 0;
    static constexpr int offset_y = 0;
    static constexpr int offset_z = 0;
};

// Base template class for block storage to reduce code duplication
template <int BlockSize, int Directions>
class BlockBase {
   public:
    BlockBase() {
    }

    // Calculate linear index from 3D coordinates
    // i - row index (0 to BlockSize-1)
    // j - column index (0 to BlockSize-1)
    // d - direction index (0 to Directions-1)
    inline double& operator()(int i, int j, int d) noexcept {
        const int index = calculateIndex(i, j, d);
        assert(index < std::ssize(values) && "Index out of bounds");
        return values[index];
    }

    const double& operator()(int i, int j, int d) const {
        const int index = calculateIndex(i, j, d);
        assert(index < std::ssize(values) && "Index out of bounds");
        return values[index];
    }

    void resize(int m) {
        values.resize(m);
    }

    void setZero() {
        std::fill(values.begin(), values.end(), 0.0);
    }

    void clear() {
        values.clear();
    }

    int size() const {
        return values.size();
    }

    std::vector<double> values;

   private:
    // Calculates linear index from 3D coordinates
    inline int calculateIndex(int i, int j, int d) const {
        assert(i >= 0 && j >= 0 && d >= 0);
        assert(i < BlockSize && j < BlockSize && d < Directions);
        return j + BlockSize * (BlockSize * d + i);
    }
};

// Template class for block matrices to reduce code duplication
template <typename BlockType, int BlockCapacity>
class BlockMatrixBase {
   public:
    BlockMatrixBase(size_t size) : data(size) {
    }
    BlockMatrixBase() {
    }

    void resize(int num_elements) {
        data.resize(num_elements);
        non_zeros.resize(num_elements);
    }

    void reserve() {
        for (auto& block : data) {
            block.resize(BlockCapacity);
        }
    }

    void setBlocksZero() {
        RECORD_TIMER;
#pragma omp parallel for schedule(dynamic, 128)
        for (unsigned long i = 0; i < non_zeros.size(); i++) {
            if (non_zeros[i]) {
                // Only resize if necessary
                if (data[i].size() != BlockCapacity) {
                    data[i].resize(BlockCapacity);
                }
                data[i].setZero();
            } else {
                data[i].clear();
            }
        }
    }

    void setZero() {
        particlesInCell(data.size() * sizeof(data[0]));
#pragma omp parallel for schedule(dynamic, 128)
        for (auto& v : data) {
            if (v.size() != BlockCapacity) {
                v.resize(BlockCapacity);
            }
            v.setZero();
        }
    }

    void get_nonzerosCells(Array3D<int>& particlesInCell) {
        RECORD_TIMER_PARAMS(particlesInCell.capacity());
#pragma omp parallel for schedule(static, 64 * 8 * 2)
        for (int i = 0; i < particlesInCell.capacity(); i++) {
            non_zeros[i] = (non_zeros[i] || particlesInCell(i) != 0);
        }
    }

    void prepare(Array3D<int>& particlesInCell) {
        RECORD_TIMER;
        std::fill(non_zeros.begin(), non_zeros.end(), false);
        get_nonzerosCells(particlesInCell);
        setBlocksZero();
    }

    BlockType& operator[](int i) {
        assert(i >= 0 && i < std::ssize(data) && "Index out of bounds");
        return data[i];
    }

    const BlockType& operator[](int i) const {
        assert(i >= 0 && i < std::ssize(data) && "Index out of bounds");
        return data[i];
    }

    std::vector<bool> non_zeros;

   private:
    std::vector<BlockType> data;
};

// Block is a part of the Lapenta matrix for a given cell (i,j,k) - only for
// Linear shape factor Block contains 12 * 12 * 9 elements d - combined
// directions - XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ i - grid index for G. Storage
// in local dense format - 3x 2y 2z for X, 2x 3y 2z for Y, 2x 2y 3z for Z j -
// grid index for G'. Storage in local dense format - 3x 2y 2z for X, 2x 3y 2z
// for Y, 2x 2y 3z for Z
using Block = BlockBase<BlockDims::BLOCK_SIZE, BlockDims::DIRECTIONS>;

// NGP version uses smaller blocks (4x4) with the same number of directions
using BlockNGP = BlockBase<BlockDimsNGP::BLOCK_SIZE, BlockDimsNGP::DIRECTIONS>;

// Define specialized block matrix types
using BlockMatrix = BlockMatrixBase<Block, BlockDims::BLOCK_CAPACITY>;
using BlockMatrixNGP = BlockMatrixBase<BlockNGP, BlockDimsNGP::BLOCK_CAPACITY>;
