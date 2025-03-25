#pragma once

#include <assert.h>

#include <unordered_map>
#include <vector>

#include "util.h"
#include "containers.h"

typedef std::unordered_map<int, double> IndexMap;

// sizeX = 3; sizeY = 2; sizeZ = 2
// inline constexpr int indX(int x, int y, int z)  {
//     return z + 2 * (y + 2 * x);
// }
// // sizeX = 2; sizeY = 3; sizeZ = 2;

// inline constexpr int indY(int x, int y, int z)  {
//     return z + 2 * (y + 3 * x);
// }
// // sizeX = 2; sizeY = 2; sizeZ = 3;

// inline constexpr int indZ(int x, int y, int z)  {
//     return z + 3 * (y + 2 * x);
// }
// inline constexpr int ind(int x, int y, int z) {
//     return z + Nz * (y + Ny * x + Nx * 0);
// }

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

constexpr int BLOCK_CAPACITY = 12 * 12 * 9;

// Block is a part of the Lapenta matrix for a given cell (i,j,k) - only for Linaer shape factor
// Block contains 12 * 12 * 9 elements
// d - combined directions - XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ
// i - grid index for G. Storage in local dense format - 3x 2y 2z for X, 2x 3y 2z for Y, 2x 2y 3z for Z
// j - grid index for G'. Storage in local dense format - 3x 2y 2z for X, 2x 3y 2z for Y, 2x 2y 3z for Z
class Block {
   public:
   Block(){}
    inline double& operator()(int i, int j, int d) noexcept {
        // return values[d + 9 * (12 * j + i)];
        // assert j + 12 * (12 * d + i) < BLOCK_CAPACITY
        // assert values.size() == BLOCK_CAPACITY
        return values[j + 12 * (12 * d + i)];
    }
    const double& operator()(int i, int j, int d) const {
        // return values[d + 9 * (12 * j + i)];
        // assert j + 12 * (12 * d + i) < BLOCK_CAPACITY
        // assert values.size() == BLOCK_CAPACITY
        return values[j + 12 * (12 * d + i)];
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
    int size() { return values.size(); }
    std::vector<double> values;
};

namespace BlockDims {
namespace X {
constexpr int I = 3, J = 2, K = 2;
}
namespace Y {
constexpr int I = 2, J = 3, K = 2;
}
namespace Z {
constexpr int I = 2, J = 2, K = 3;
}

constexpr int BLOCK = 12;
constexpr int DIRS = 9;
constexpr int ELEMS = BLOCK * BLOCK * DIRS;
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

// Структура для хранения элемента разреженной матрицы
class Triplet {

public:
    Triplet(int r, int c, double v) : _row(r), _col(c), _value(v) {}
    Triplet() : _row(0), _col(0), _value(0) {}
    const int& row() const noexcept { return _row; }
    const int& col() const noexcept { return _col; }
    const double& value() const noexcept { return _value; }
    int& row() {return _row;}
    int& col() {return _col;}
    double& value() {return _value;}
    // Triplet(Triplet&& other) noexcept = default;
    //Triplet& operator=(Triplet&& other) noexcept = default;

    private:
     int _row;
     int _col;
     double _value;
};

template <typename RowIdx, typename ColIdx, int DIR>
void processComponent(int i_cell, int j_cell, int k_cell, const Block& block,
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
// Функция, которая заполняет матрицу напрямую через внутренние массивы Eigen.
// Предполагается, что:
//   - trips отсортирован по (row, col)
//   - вектор trips не содержит дубликатов
//   - матрица хранится в формате RowMajor
static void optimizedSetFromTriplets(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                              const std::vector<Triplet>& trips) {
    int numRows = mat.rows();
                              int numCols = mat.cols();
    const int nnz = trips.size();

    // Устанавливаем размеры матрицы и выделяем ровно nnz элементов
   // mat.resize(numRows, numCols);
    mat.resizeNonZeros(nnz);

    // Получаем указатели на внутренние массивы
    int* outer = mat.outerIndexPtr();
    int* inner = mat.innerIndexPtr();
    double* values = mat.valuePtr();

    // 1. Подсчитываем число ненулевых элементов в каждой строке
    std::vector<int> rowCounts(numRows, 0);
    // Также вычисляем для каждой строки индекс первого элемента в trips.
    // Поскольку trips отсортирован по строкам, эти данные можно вычислить одним
    // проходом.
    std::vector<int> rowStart(numRows, -1);
    for (int i = 0; i < nnz; ++i) {
        int r = trips[i].row();
        ++rowCounts[r];
        if (rowStart[r] == -1)
            rowStart[r] = i;
    }

    // 2. Строим массив outer как префиксную сумму rowCounts
    outer[0] = 0;
    for (int i = 0; i < numRows; ++i) {
        outer[i + 1] = outer[i] + rowCounts[i];
    }

// 3. Заполняем внутренние массивы inner и values
// Распараллеливаем цикл по строкам, так как для каждой строки область
// заполнения определяется outer.
#pragma omp parallel for schedule(static)
    for (int r = 0; r < numRows; ++r) {
        if (rowCounts[r] == 0)
            continue;
        int start = rowStart[r];   // индекс первого триплета для строки r
        for (int j = 0; j < rowCounts[r]; ++j) {
            inner[outer[r] + j] = trips[start + j].col();
            values[outer[r] + j] = trips[start + j].value();
        }
    }

    // 4. Завершаем заполнение матрицы: объявляем, что данные уже сжаты.
    mat.makeCompressed();
}

struct pair_hash {
    size_t operator()(const std::pair<int, int>& p) const {
        return static_cast<size_t>(p.first) << 32 | p.second;
    }
};

typedef std::vector<Block> BMatrix2;

class BlockMatrix {
   public:

    BlockMatrix(size_t size) : data(size) {}
    BlockMatrix() {}
    void resize(int num_elements) {
        data.resize(num_elements);
        non_zeros.resize(num_elements);
    }
    void reserve(){
        for (auto& block : data) {
            block.resize(BLOCK_CAPACITY);
        }
    }

    void setBlocksZero() {
#pragma omp parallel for schedule(dynamic, 128)
        for (unsigned long i = 0; i < non_zeros.size(); i++) {
            if (non_zeros[i]) {
                data[i].resize(BLOCK_CAPACITY);
                data[i].setZero();
            } else{
                data[i].clear();
            }

        }
    }
        void setZero() {
#pragma omp parallel for schedule(dynamic, 128)
        for (auto& v : data) {
                v.resize(BLOCK_CAPACITY);
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

    Block& operator[](int i) { return data[i]; }
    const Block& operator[](int i) const { return data[i]; }
    std::vector<bool> non_zeros;

   private:
    std::vector<Block> data;
};

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
