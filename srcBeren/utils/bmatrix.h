#pragma once

#include <assert.h>

#include <unordered_map>
#include <vector>

#include "util.h"
#include "containers.h"
#include "service.h"

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

// Common constants
constexpr int BLOCK_SIZE = 12;
constexpr int DIRECTIONS = 9;
constexpr int BLOCK_CAPACITY = BLOCK_SIZE * BLOCK_SIZE * DIRECTIONS;
}   // namespace BlockDims

// Improved indexing functions with named constants
inline constexpr int indX(int i, int j, int k) {
    return i * (BlockDims::X::SIZE_J * BlockDims::X::SIZE_K) +
           j * BlockDims::X::SIZE_K + k;
}

inline constexpr int indY(int i, int j, int k) {
    return i * (BlockDims::Y::SIZE_J * BlockDims::Y::SIZE_K) +
           j * BlockDims::Y::SIZE_K + k;
}

inline constexpr int indZ(int i, int j, int k) {
    return i * (BlockDims::Z::SIZE_J * BlockDims::Z::SIZE_K) +
           j * BlockDims::Z::SIZE_K + k;
}

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
    return i * (X::SIZE_J * X::SIZE_K) +
           j * X::SIZE_K + k;
}

inline constexpr int indY(int i, int j, int k) {
    return i * (Y::SIZE_J * Y::SIZE_K) +
           j * Y::SIZE_K + k;
}

inline constexpr int indZ(int i, int j, int k) {
    return i * (Z::SIZE_J * Z::SIZE_K) +
           j * Z::SIZE_K + k;
}

// Common constants
constexpr int BLOCK_SIZE = 4;
constexpr int DIRECTIONS = 9;
constexpr int BLOCK_CAPACITY = BLOCK_SIZE * BLOCK_SIZE * DIRECTIONS;
}   // namespace BlockDims



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
template<int BlockSize, int Directions>
class BlockBase {
public:
    BlockBase() {}
    
    // Calculate linear index from 3D coordinates
    // i - row index (0 to BlockSize-1)
    // j - column index (0 to BlockSize-1)
    // d - direction index (0 to Directions-1)
    inline double& operator()(int i, int j, int d) noexcept {
        const int index = calculateIndex(i, j, d);
        assert(index < values.size() && "Index out of bounds");
        return values[index];
    }
    
    const double& operator()(int i, int j, int d) const {
        const int index = calculateIndex(i, j, d);
        assert(index < values.size() && "Index out of bounds");
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
        return j + BlockSize * (BlockSize * d + i);
    }
};

// Block is a part of the Lapenta matrix for a given cell (i,j,k) - only for Linear shape factor
// Block contains 12 * 12 * 9 elements
// d - combined directions - XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ
// i - grid index for G. Storage in local dense format - 3x 2y 2z for X, 2x 3y 2z for Y, 2x 2y 3z for Z
// j - grid index for G'. Storage in local dense format - 3x 2y 2z for X, 2x 3y 2z for Y, 2x 2y 3z for Z
using Block = BlockBase<BlockDims::BLOCK_SIZE, BlockDims::DIRECTIONS>;

// NGP version uses smaller blocks (4x4) with the same number of directions
using BlockNGP = BlockBase<BlockDimsNGP::BLOCK_SIZE, BlockDimsNGP::DIRECTIONS>;

// Template class for block matrices to reduce code duplication
template<typename BlockType, int BlockCapacity>
class BlockMatrixBase {
public:
    BlockMatrixBase(size_t size) : data(size) {}
    BlockMatrixBase() {}
    
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
        #pragma omp parallel for schedule(dynamic, 128)
        for (auto& v : data) {
            if (v.size() != BlockCapacity) {
                v.resize(BlockCapacity);
            }
            v.setZero();
        }
    }

    void get_nonzerosCells(Array3D<int>& particlesInCell) {
        for (int i = 0; i < particlesInCell.capacity(); i++) {
            non_zeros[i] = (non_zeros[i] || particlesInCell(i) != 0);
        }
    }
    
    void prepare(Array3D<int>& particlesInCell) {
        std::fill(non_zeros.begin(), non_zeros.end(), false);
        get_nonzerosCells(particlesInCell);
        setBlocksZero();
    }

    BlockType& operator[](int i) { 
        assert(i >= 0 && i < data.size() && "Index out of bounds");
        return data[i]; 
    }
    
    const BlockType& operator[](int i) const { 
        assert(i >= 0 && i < data.size() && "Index out of bounds");
        return data[i]; 
    }
    
    std::vector<bool> non_zeros;

private:
    std::vector<BlockType> data;
};

// Define specialized block matrix types
using BlockMatrix = BlockMatrixBase<Block, BlockDims::BLOCK_CAPACITY>;
using BlockMatrixNGP = BlockMatrixBase<BlockNGP, BlockDimsNGP::BLOCK_CAPACITY>;

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
void processComponent(int i_cell, int j_cell, int k_cell, const Block_t& block,
                      std::vector<Triplet>& trips, [[maybe_unused]] int Nx,
                      int Ny, int Nz, double tolerance) {
    auto vind = [&](int i, int j, int k, int d) {
        return d + 3 * (i * Ny * Nz + j * Nz + k);
    };
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1)
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1)
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row =
                    vind(i_cell + x1 + RowIdx::offset_x,
                         j_cell + y1 + RowIdx::offset_y,
                         k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                for (int x2 = 0; x2 < ColIdx::size_x; ++x2)
                    for (int y2 = 0; y2 < ColIdx::size_y; ++y2)
                        for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                            const double val =
                                block(RowIdx::calculate(x1, y1, z1),
                                      ColIdx::calculate(x2, y2, z2), DIR);

                            if (std::abs(val) > tolerance) {
                                const int col =
                                    vind(i_cell + x2 + ColIdx::offset_x,
                                         j_cell + y2 + ColIdx::offset_y,
                                         k_cell + z2 + ColIdx::offset_z,
                                         ColIdx::dir);

                                trips.emplace_back(Triplet{row, col, val});
                            }
                        }
            }
}

typedef std::vector<Block> BMatrix2;

class MatrixCSR {
   private:
    std::vector<double> values;
    std::vector<int> col_indices;
    std::vector<int> row_ptrs;
    int num_rows;
    int num_cols;

   public:
    MatrixCSR(int rows, int cols) : num_rows(rows), num_cols(cols) {
        row_ptrs.resize(rows + 1, 0);
    }

    void reserve(size_t nnz) {
        values.reserve(nnz);
        col_indices.reserve(nnz);
    }

    const std::vector<double>& get_values() const { return values; }
    const std::vector<int>& get_col_indices() const { return col_indices; }
    const std::vector<int>& get_row_ptrs() const { return row_ptrs; }
    std::vector<int>& get_row_ptrs() { return row_ptrs; }

    std::vector<double> multiply(const std::vector<double>& x) const;

};
