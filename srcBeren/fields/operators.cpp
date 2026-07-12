#include "operators.h"

// #include <execution>
#include <atomic>

#include "Mesh.h"
#include "Shape.h"
#include "World.h"
#include "config.h"
#include "log_macros.h"
#include "pmms.hpp"
#include "timer.h"
#include "util.h"

// matrix ColMajor
// (0,0) (0,1) (0,2)
// (1,0) (1,1) (1,2)
// (2,0) (2,1) (2,2)

// Ex(i+/2,j,k), Ey(i,j+1/2,k), Ez(i,j,k+1/2)
// Bx(i,j+1/2,k+1/2), By(i+1/2,j,k+1/2), Bz(i+/2,j+1/2,k)

void Mesh::stencil_Imat(Operator& mat, const Domain& domain) {
    // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z();
    trips.reserve(totalSize);

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 0), vind(i, j, k, 0), 1.0));

                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 1), vind(i, j, k, 1), 1.0));

                // i,j,k
                trips.push_back(Trip(vind(i, j, k, 2), vind(i, j, k, 2), 1.0));
            }
        }
    }
    LOG_STEP("Stencil_Imat done. size: " << size << ". Elements: " << 3 * totalSize << " " << Imat.rows() << " "
                                         << Imat.cols() << "\n");
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_Lmat(Operator& mat, const Domain& domain) {
    //     std::vector<Trip> trips;
    //     const auto size = domain.size();
    //     const int rowsCount = 3 * size.x() * size.y() * size.z();
    //     const size_t totalSize = (size_t)rowsCount *
    //     LMAT_MAX_ELEMENTS_PER_ROW; std::cout << totalSize << "\n";
    //     trips.reserve(totalSize);

    // #pragma omp parallel
    //     {
    //         std::vector<Trip> local_trips;
    //         local_trips.reserve(totalSize / omp_get_max_threads());

    // #pragma omp for schedule(dynamic, 8)
    //         for (int row = 0; row < rowsCount; row++) {
    //             for (const auto &[col, value] : LmatX[row]) {
    //                 if (std::abs(value) > LMAT_VALUE_TOLERANCE)
    //                     local_trips.emplace_back(row, col, value);
    //             }
    //         }

    // #pragma omp critical
    //         trips.insert(trips.end(), local_trips.begin(),
    //         local_trips.end());
    //     }
    //     std::cout << "trips size: " << trips.size() << std::endl;

    //     mat.setFromTriplets(trips.begin(), trips.end());
}

static bool equalVecsTriplets(const std::vector<Triplet>& a, const std::vector<Triplet>& b) {
    if (a.size() != b.size()) {
        std::cout << "size a " << a.size() << " size b " << b.size() << "\n";
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].col() != b[i].col() || a[i].row() != b[i].row() || fabs(a[i].value() - b[i].value()) > 1.e-15) {
            std::cout << "row " << a[i].row() << " col " << a[i].col() << " value " << a[i].value() << " row "
                      << b[i].row() << " col " << b[i].col() << " value " << b[i].value() << " "
                      << std::abs(a[i].value() - b[i].value()) << "\n";
            return false;
        }
    return true;
}

std::vector<Triplet> multyPhaseMerge(std::vector<std::vector<Triplet>>& local_vectors) {
    RECORD_TIMER;

    // Собираем только непустые локальные векторы для дальнейшего слияния
    std::vector<std::vector<Triplet>> non_empty;
    for (auto& v : local_vectors) {
        if (!v.empty()) {
            non_empty.push_back(std::move(v));
        }
    }
    if (non_empty.empty()) {
        return {};
    }

    // Многофазное слияние: объединяем пары отсортированных векторов параллельно
    // с предварительным резервированием памяти
    while (non_empty.size() > 1) {
        size_t new_size = (non_empty.size() + 1) / 2;
        std::vector<std::vector<Triplet>> new_vectors(new_size);
        size_t pairs = non_empty.size() / 2;

#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < pairs; ++i) {
            timer::flatTimer timerIt("single merge");
            const auto& left = non_empty[2 * i];
            const auto& right = non_empty[2 * i + 1];
            // Резервируем память для слияния двух векторов
            size_t merged_capacity = left.size() + right.size();
            std::vector<Triplet> merged;
            merged.reserve(merged_capacity);

            size_t li = 0, ri = 0;
            while (li < left.size() && ri < right.size()) {
                if (compareTriplets(left[li], right[ri])) {
                    merged.push_back(left[li]);
                    ++li;
                } else if (compareTriplets(right[ri], left[li])) {
                    merged.push_back(right[ri]);
                    ++ri;
                } else {
                    // Если ключи равны – складываем значения
                    Triplet t = left[li];
                    t.value() += right[ri].value();
                    merged.push_back(t);
                    ++li;
                    ++ri;
                }
            }
            while (li < left.size()) {
                merged.push_back(left[li++]);
            }
            while (ri < right.size()) {
                merged.push_back(right[ri++]);
            }
            new_vectors[i] = std::move(merged);
        }
        // Если число векторов нечётное – последний переносим без изменений
        if (non_empty.size() % 2 == 1) {
            new_vectors.back() = std::move(non_empty.back());
        }
        non_empty = std::move(new_vectors);
    }
    return std::move(non_empty[0]);
}

template <int maxNnz = 12 * 12 * 9>
struct RowInfo {
    RowInfo(int rowIn) : row(rowIn), nnz(0) {
    }
    RowInfo() {
    }

    ~RowInfo() {
    }

    template <int otherNnz>
    RowInfo(int count, const RowInfo<otherNnz>* others) {
        mergeFromOthers(count, others);
    }

    RowInfo& operator=(const RowInfo<maxNnz>& other) {
        nnz = other.nnz;
        row = other.row;
        assert(nnz <= maxNnz);

        std::copy_n(other.values.begin(), nnz, values.begin());
        std::copy_n(other.columns.begin(), nnz, columns.begin());

        return *this;
    }

    RowInfo(const RowInfo& other) {
        nnz = other.nnz;
        row = other.row;

        std::copy_n(other.values.begin(), nnz, values.begin());
        std::copy_n(other.columns.begin(), nnz, columns.begin());
    }

    bool operator!=(const RowInfo& other) {
        if (row != other.row || nnz != other.nnz) {
            return true;
        }
        for (int i = 0; i < nnz; ++i) {
            if (values[i] != other.values[i] || columns[i] != other.columns[i]) {
                return true;
            }
        }
        return false;
    }

    void addValue(int col, double val) {
        for (int i = 0; i < nnz; ++i) {
            if (columns[i] == col) {
                values[i] += val;
                return;
            }
        }
        assert(nnz < maxNnz);
        values[nnz] = val;
        columns[nnz] = col;
        nnz += 1;
    }

    void sort() {
        std::array<std::pair<int, double>, maxNnz> pairs;
        for (int i = 0; i < nnz; ++i) {
            pairs[i] = {columns[i], values[i]};
        }

        std::sort(pairs.begin(), pairs.begin() + nnz,
                  [](const std::pair<int, double>& a, const std::pair<int, double>& b) { return a.first < b.first; });

        for (int i = 0; i < nnz; ++i) {
            columns[i] = pairs[i].first;
            values[i] = pairs[i].second;
        }
    }

    template <int otherNnz>
    void mergeFromOthers(int count, const RowInfo<otherNnz>* others) {
        // static int maxCount = count;
        // maxCount = std::max(count, maxCount);
        // std::cout<<"max count: "<<maxCount<<std::endl;

        int its[count];
        // more cache-friendly access
        int othersNnz[count];
        std::fill_n(its, count, 0);

        for (int i = 0; i < count; ++i) {
            othersNnz[i] = others[i].nnz;
        }

        for (int i = 1; i < count; ++i) {
            assert(others[i].row == others[0].row);
        }

        row = others[0].row;
        nnz = 0;

        while (true) {
            int smallestCol = std::numeric_limits<int>::max();
            for (int i = 0; i < count; ++i) {
                if (its[i] < othersNnz[i])
                    smallestCol = std::min(smallestCol, others[i].columns[its[i]]);
            }

            if (smallestCol == std::numeric_limits<int>::max()) {
                break;
            }

            double acc = 0.0;
            for (int i = 0; i < count; ++i) {
                if (its[i] < othersNnz[i] && others[i].columns[its[i]] == smallestCol) {
                    acc += others[i].values[its[i]];
                    its[i] += 1;
                }
            }

            assert(nnz < maxNnz);
            columns[nnz] = smallestCol;
            values[nnz] = acc;
            nnz += 1;
        }
    }

    int row;
    int nnz;
    std::array<double, maxNnz> values;
    std::array<int, maxNnz> columns;
};

struct RowInfosCompressed {
    RowInfosCompressed() {
        outer.reserve(1024);
        values.reserve(1024);
        cols.reserve(1024);

        outer.push_back(0);
    }

    template <int maxNnz>
    void mergeFromRowInfo(int count, const RowInfo<maxNnz>* rowInfos) {
        for (int i = 1; i < count; ++i) {
            assert(rowInfos[i].row == rowInfos[0].row);
        }

        const int row = rowInfos[0].row;
        for (; prevRow < row; ++prevRow) {
            outer.push_back(outer.back());
        }

        int its[count];
        std::fill_n(its, count, 0);
        while (true) {
            int smallestCol = std::numeric_limits<int>::max();
            for (int i = 0; i < count; ++i) {
                if (its[i] < rowInfos[i].nnz)
                    smallestCol = std::min(smallestCol, rowInfos[i].columns[its[i]]);
            }

            if (smallestCol == std::numeric_limits<int>::max()) {
                break;
            }

            double acc = 0.0;
            for (int i = 0; i < count; ++i) {
                if (its[i] < rowInfos[i].nnz && rowInfos[i].columns[its[i]] == smallestCol) {
                    acc += rowInfos[i].values[its[i]];
                    its[i] += 1;
                }
            }

            cols.push_back(smallestCol);
            values.push_back(acc);
            outer.back() += 1;
        }
    }

    std::vector<int64_t> outer;
    std::vector<double> values;
    std::vector<int64_t> cols;
    int prevRow = 0;
};

template <typename ColIdx, int DIR, typename Block_t>
void processRow(int i_cell, int j_cell, int k_cell, const Block_t& block, int rowInBlock, [[maybe_unused]] int Nx,
                int Ny, int Nz, double tolerance, RowInfo<>& rowInfo) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };

    for (int x2 = 0; x2 < ColIdx::size_x; ++x2)
        for (int y2 = 0; y2 < ColIdx::size_y; ++y2)
            for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                const double val = block(rowInBlock, ColIdx::calculate(x2, y2, z2), DIR);

                if (std::abs(val) > tolerance) {
                    const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                         k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                    rowInfo.addValue(col, val);
                }
            }
}
template <typename RowIdx, int zeroDir>
void processRow2(const Mesh& mesh, int row, int max_i, int max_j, int max_k, int border, [[maybe_unused]] int xSize,
                 int ySize, int zSize, double tolerance, RowInfo<>& rowInfo) {
    const int k0 = (row / 3) % zSize - RowIdx::offset_z;
    const int j0 = (row / 3 / zSize) % ySize - RowIdx::offset_y;
    const int i0 = (row / 3 / zSize / ySize) - RowIdx::offset_x;

    // X component
    for (int i = std::max(border, i0 - RowIdx::size_x + 1); i <= std::min(i0, max_i); ++i) {
        for (int j = std::max(border, j0 - RowIdx::size_y + 1); j <= std::min(j0, max_j); ++j) {
            for (int k = std::max(border, k0 - RowIdx::size_z + 1); k <= std::min(k0, max_k); ++k) {
                if (!mesh.LmatX2.non_zeros[mesh.sind(i, j, k)])
                    continue;

                const Block& block = mesh.LmatX2[mesh.sind(i, j, k)];

                const int rowInBlock = RowIdx::calculate(i0 - i, j0 - j, k0 - k);

                processRow<XIndexer, zeroDir + 0>(i, j, k, block, rowInBlock, xSize, ySize, zSize, tolerance, rowInfo);
                processRow<YIndexer, zeroDir + 1>(i, j, k, block, rowInBlock, xSize, ySize, zSize, tolerance, rowInfo);
                processRow<ZIndexer, zeroDir + 2>(i, j, k, block, rowInBlock, xSize, ySize, zSize, tolerance, rowInfo);
            }
        }
    }

    // const int k0 = (row / 3) % zSize;
    // const int j0 = (row / 3 / zSize) % ySize;
    // const int i0 = (row / 3 / zSize / ySize);

    // constexpr int commonSize = 2;

    // // X component
    // for (int i = std::max(border, i0 - commonSize + 1); i <= std::min(i0, max_i); ++i) {
    //     for (int j = std::max(border, j0 - commonSize + 1); j <= std::min(j0, max_j); ++j) {
    //         for (int k = std::max(border, k0 - commonSize + 1); k <= std::min(k0, max_k); ++k) {
    //             if (!mesh.LmatX2.non_zeros[mesh.sind(i, j, k)])
    //                 continue;

    //             const Block& block = mesh.LmatX2[mesh.sind(i, j, k)];

    //             // TODO: add other substraction
    //             const int rowInBlockX = RowIdx::calculate(i0 - i - XIndexer::offset_x, j0 - j, k0 - k);
    //             const int rowInBlockY = RowIdx::calculate(i0 - i, j0 - j - YIndexer::offset_y, k0 - k);
    //             const int rowInBlockZ = RowIdx::calculate(i0 - i, j0 - j, k0 - k - ZIndexer::offset_z);

    //             processRow<XIndexer, 0>(i, j, k, block, rowInBlockX, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<YIndexer, 1>(i, j, k, block, rowInBlockX, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<ZIndexer, 2>(i, j, k, block, rowInBlockX, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<XIndexer, 3>(i, j, k, block, rowInBlockY, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<YIndexer, 4>(i, j, k, block, rowInBlockY, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<ZIndexer, 5>(i, j, k, block, rowInBlockY, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<XIndexer, 6>(i, j, k, block, rowInBlockZ, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<YIndexer, 7>(i, j, k, block, rowInBlockZ, xSize, ySize, zSize, tolerance, rowInfo);
    //             processRow<ZIndexer, 8>(i, j, k, block, rowInBlockZ, xSize, ySize, zSize, tolerance, rowInfo);
    //         }
    //     }
    // }
}

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
void processComponent2(int i_cell, int j_cell, int k_cell, const Block_t& block, [[maybe_unused]] int Nx, int Ny,
                       int Nz, double tolerance, std::vector<RowInfo<12>>& rowInfos) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row = vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                     k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                bool isAddedInThisRow = false;

                const int rowIdx = RowIdx::calculate(x1, y1, z1);

                for (int x2 = 0; x2 < ColIdx::size_x; ++x2) {
                    for (int y2 = 0; y2 < ColIdx::size_y; ++y2) {
                        for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                            const int colIdx = ColIdx::calculate(x2, y2, z2);
                            const double val = block(rowIdx, colIdx, DIR);

                            // std::cout << "test coords and val: " << x1 << " " << y1 << " " << z1 << " " << x2 << " "
                            //           << y2 << " " << z2 << " " << val << " for row " << row << std::endl;

                            if (std::abs(val) > tolerance) {
                                // std::cout << "test add to row " << row << std::endl;
                                const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                                     k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                                if (!isAddedInThisRow) [[unlikely]] {
                                    rowInfos.emplace_back(row);
                                    isAddedInThisRow = true;
                                }
                                rowInfos.back().addValue(col, val);
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
void computeNonZeroBlocksInRows(int i_cell, int j_cell, int k_cell, const Block_t& block, [[maybe_unused]] int Nx,
                                int Ny, int Nz, double tolerance, std::vector<uint8_t>& nonZeroBlocks) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row = vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                     k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                bool isAddedInThisRow = false;

                const int rowIdx = RowIdx::calculate(x1, y1, z1);
                for (int x2 = 0; x2 < ColIdx::size_x && !isAddedInThisRow; ++x2) {
                    for (int y2 = 0; y2 < ColIdx::size_y && !isAddedInThisRow; ++y2) {
                        for (int z2 = 0; z2 < ColIdx::size_z && !isAddedInThisRow; ++z2) {
                            const int colIdx = ColIdx::calculate(x2, y2, z2);
                            const double val = block(rowIdx, colIdx, DIR);

                            if (std::abs(val) > tolerance) {
                                std::atomic_ref<uint8_t> toUpd(nonZeroBlocks[row]);
                                ++toUpd;
                                isAddedInThisRow = true;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
void processComponent3(int i_cell, int j_cell, int k_cell, const Block_t& block, [[maybe_unused]] int Nx, int Ny,
                       int Nz, double tolerance, std::vector<RowInfo<12>>& rowInfosGlobal,
                       std::vector<uint8_t>& nonZeroBlocksCount, const std::vector<int>& blockStarts) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };

    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row = vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                     k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                bool isAddedInThisRow = false;

                const int rowIdx = RowIdx::calculate(x1, y1, z1);
                int blockIx = -1;

                for (int x2 = 0; x2 < ColIdx::size_x; ++x2) {
                    for (int y2 = 0; y2 < ColIdx::size_y; ++y2) {
                        for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                            const int colIdx = ColIdx::calculate(x2, y2, z2);
                            const double val = block(rowIdx, colIdx, DIR);

                            if (std::abs(val) > tolerance) {
                                // std::cout << "test add to row " << row << std::endl;
                                const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                                     k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                                if (!isAddedInThisRow) [[unlikely]] {
                                    isAddedInThisRow = true;
                                    std::atomic_ref<uint8_t> blockIxRef(nonZeroBlocksCount[row]);
                                    blockIx = (blockIxRef++) + blockStarts[row];
                                    rowInfosGlobal[blockIx].nnz = 0;
                                    rowInfosGlobal[blockIx].row = row;
                                }
                                rowInfosGlobal[blockIx].addValue(col, val);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Mesh::stencil_Lmat2_OPT1(Operator& mat, const Domain& domain) const {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    // static int check_count = 0;
    // std::vector<Triplet> trips;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;
    // const int num_threads = std::min(omp_get_max_threads(), 128);

    const int rows = mat.rows();
    const int nthr = omp_get_max_threads();
    // const int nthr = 1;

    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * ySize * zSize + j * zSize + k); };

    timer::commonTimer timerCountingRows("getting non zero rows");
    // std::vector<int> nonEmptyRows(1024);
    // // #pragma omp for collapse(3) schedule(dynamic, 8) nowait
    // for (int i = BORDER; i < max_i; ++i) {
    //     for (int j = BORDER; j < max_j; ++j) {
    //         for (int k = BORDER; k < max_k; ++k) {
    //             if (LmatX2.non_zeros[sind(i, j, k)]) {
    //                 for (int i2 = -1; i2 < 2; ++i2) {
    //                     for (int j2 = -1; j2 < 2; ++j2) {
    //                         for (int k2 = -1; k2 < 2; ++k2) {
    //                             nonEmptyRows.push_back(vind(i + i2, j + j2, k + k2, 0));
    //                             nonEmptyRows.push_back(vind(i + i2, j + j2, k + k2, 1));
    //                             nonEmptyRows.push_back(vind(i + i2, j + j2, k + k2, 2));
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // timerCountingRows.finish();

    std::vector<std::vector<RowInfo<>>> localRowInfos;
    localRowInfos.resize(nthr);

    for (int i = 0; i < nthr; ++i) {
        localRowInfos[i].reserve(1024);
        localRowInfos[i].emplace_back(-1);
    }

    timer::commonTimer timerUnpacking("unpacking");
    // #pragma omp parallel for num_threads(nthr)
    for (int row = 0; row < rows; ++row) {
        // for (int i = 0; i < std::ssize(nonEmptyRows); ++i) {
        // const int row = nonEmptyRows[i];

        std::vector<RowInfo<>>& rowInfos = localRowInfos[omp_get_thread_num()];

        if (rowInfos.back().nnz != 0) {
            rowInfos.emplace_back(row);
        }

        RowInfo<>& rowInfo = rowInfos.back();
        rowInfo.row = row;

        // X component
        if (row % 3 == 0)
            processRow2<XIndexer, 0>(*this, row, max_i, max_j, max_k, BORDER, xSize, ySize, zSize, TOL, rowInfo);
        // Y component
        if (row % 3 == 1)
            processRow2<YIndexer, 3>(*this, row, max_i, max_j, max_k, BORDER, xSize, ySize, zSize, TOL, rowInfo);
        // Z component
        if (row % 3 == 2)
            processRow2<ZIndexer, 6>(*this, row, max_i, max_j, max_k, BORDER, xSize, ySize, zSize, TOL, rowInfo);

        if (rowInfo.nnz != 0) {
            rowInfo.sort();
        }
    }
    timerUnpacking.finish();

    int totalNnz = 0;

    std::vector<int> outerIndexes(rows + 1);
    std::fill(outerIndexes.begin(), outerIndexes.end(), 0);

#pragma omp parallel for reduction(+ : totalNnz)
    for (int i = 0; i < nthr; ++i) {
        const std::vector<RowInfo<>>& rowInfos = localRowInfos[i];
        for (int j = 0; j < std::ssize(rowInfos); ++j) {
            outerIndexes[rowInfos[j].row + 1] = rowInfos[j].nnz;
            totalNnz += rowInfos[j].nnz;
        }
    }

    for (int i = 1; i < rows + 1; ++i) {
        outerIndexes[i] += outerIndexes[i - 1];
    }

    timer::commonTimer timerResize("mat.resizeNonZeros()");
    mat.resizeNonZeros(totalNnz);
    timerResize.finish();

    int* outer = mat.outerIndexPtr();
    int* ind = mat.innerIndexPtr();
    // int* innerNnzRes = res.innerNonZeroPtr();

    double* values = mat.valuePtr();

    timer::commonTimer timerFilling("filling");

    outer[0] = 0;
    int ix1 = 0;
    int ix2 = 0;

    while (localRowInfos[ix1].size() == 0) {
        ix1 += 1;
    }

    for (int i = 0; i < rows; ++i) {
        outer[i + 1] = outerIndexes[i + 1];

        const int start = outerIndexes[i];
        const int end = outerIndexes[i + 1];
        if (start == end)
            continue;

        const RowInfo rowInfo = localRowInfos[ix1][ix2];
        assert(rowInfo.row == i);
        for (int j = start; j < end; ++j) {
            values[j] = rowInfo.values[j - start];
            ind[j] = rowInfo.columns[j - start];
        }
        ix2 += 1;
        if (ix2 == std::ssize(localRowInfos[ix1])) {
            ix1 += 1;
            ix2 = 0;
        }
    }
}

void Mesh::stencil_Lmat2_OPT2(Operator& mat, const Domain& domain) const {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    // static int check_count = 0;
    // std::vector<Triplet> trips;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;
    // const int num_threads = std::min(omp_get_max_threads(), 128);

    const int rows = mat.rows();
    const int nthr = omp_get_max_threads();
    // const int nthr = 1;

    std::vector<std::vector<RowInfo<12>>> rowInfosTest(nthr);

    timer::commonTimer timerUnpackingNew("unpacking new");
#pragma omp parallel num_threads(nthr)
    {
        std::vector<RowInfo<12>>& localRowInfos = rowInfosTest[omp_get_thread_num()];
        localRowInfos.reserve(1024 * 1024 * 4);
        timer::commonTimer timerOmp("OMP loop");
#pragma omp for collapse(3) schedule(dynamic, 8)
        for (int i = BORDER; i < max_i; ++i) {
            for (int j = BORDER; j < max_j; ++j) {
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = LmatX2[sind(i, j, k)];

                    // X component
                    processComponent2<XIndexer, XIndexer, 0>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<XIndexer, YIndexer, 1>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<XIndexer, ZIndexer, 2>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);

                    // Y component
                    processComponent2<YIndexer, XIndexer, 3>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<YIndexer, YIndexer, 4>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<YIndexer, ZIndexer, 5>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);

                    // Z component
                    processComponent2<ZIndexer, XIndexer, 6>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<ZIndexer, YIndexer, 7>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<ZIndexer, ZIndexer, 8>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                }
            }
        }
    }
    timerUnpackingNew.finish();

    timer::commonTimer timerSortNew("sort new");
    // Indexes of each arrays is: row, orig thread owner and orig. number in its buffer

    int totalUnsortedRows = 0;
    for (int i = 0; i < nthr; ++i) {
        totalUnsortedRows += std::ssize(rowInfosTest[i]);
    }

    std::vector<int> offsets(nthr);
    offsets[0] = 0;
    for (int i = 1; i < nthr; ++i) {
        offsets[i] = std::ssize(rowInfosTest[i - 1]) + offsets[i - 1];
    }

    std::vector<std::array<int, 3>> rowPositions(totalUnsortedRows);

#pragma omp parallel for num_threads(nthr)
    for (int i = 0; i < nthr; ++i) {
        const std::vector<RowInfo<12>>& localRowInfos = rowInfosTest[i];

        for (int j = 0; j < std::ssize(localRowInfos); ++j) {
            rowPositions[offsets[i] + j] = {localRowInfos[j].row, i, j};
        }
    }

    timer::commonTimer timerStdSort("std::sort", rowPositions.size());
    std::sort(rowPositions.begin(), rowPositions.end(),
              [](const std::array<int, 3>& a, const std::array<int, 3>& b) { return a[0] < b[0]; });
    timerStdSort.finish();

    timer::commonTimer timerAllocation("allocation globalRowInfos");
    std::vector<RowInfo<12>> globalRowInfos(rowPositions.size());
    timerAllocation.finish();
    // timer::commonTimer timerAllocation2("allocation 2");
    // volatile std::vector<double> tmpVector(sizeof(RowInfo<12 * 9>));
    // timerAllocation2.finish();

#pragma omp parallel
    {
        timer::commonTimer timerOmp("OMP loop");
#pragma omp for
        for (int i = 0; i < std::ssize(rowPositions); ++i) {
            const auto [row, thread, index] = rowPositions[i];
            globalRowInfos[i] = rowInfosTest[thread][index];
        }
        timerOmp.finish();
    }

    timerSortNew.finish();

    timer::commonTimer timerMerge("merge");

    std::vector<int> blocksStarts;
    blocksStarts.reserve(std::ssize(globalRowInfos));
    blocksStarts.push_back(0);

    for (int i = 0, j = 0; i != std::ssize(globalRowInfos); i = j) {
        while (j < std::ssize(globalRowInfos) && globalRowInfos[i].row == globalRowInfos[j].row) {
            j += 1;
        }
        if (j != i) {
            blocksStarts.push_back(j);
        }
    }

    std::vector<RowInfo<12 * 12>> globalRowInfosMerged(blocksStarts.size() - 1);

#pragma omp parallel for
    for (int i = 0; i < std::ssize(blocksStarts) - 1; ++i) {
        globalRowInfosMerged[i].mergeFromOthers(blocksStarts[i + 1] - blocksStarts[i],
                                                &globalRowInfos[blocksStarts[i]]);
    }
    timerMerge.finish();

    // {
    //     timer::commonTimer testMergeTimer("test merge");

    //     RowInfosCompressed rowInfosCompressed;
    //     for (int i = 0, j = 0; i != std::ssize(globalRowInfos); i = j) {
    //         while (globalRowInfos[i].row == globalRowInfos[j].row) {
    //             j += 1;
    //         }

    //         assert(j - i <= 12 * 9 && j - i > 0);
    //         rowInfosCompressed.mergeFromRowInfo(j - i, &globalRowInfos[i]);
    //     }

    //     testMergeTimer.finish();
    // }

    timer::commonTimer timerAux("aux");
    int totalNnz = 0;
    std::vector<int> outerIndexes(rows + 1);
#pragma omp parallel for
    for (int i = 0; i < rows + 1; ++i) {
        outerIndexes[i] = 0;
    }

#pragma omp parallel for reduction(+ : totalNnz)
    for (int j = 0; j < std::ssize(globalRowInfosMerged); ++j) {
        outerIndexes[globalRowInfosMerged[j].row + 1] = globalRowInfosMerged[j].nnz;
        totalNnz += globalRowInfosMerged[j].nnz;
    }

    timer::commonTimer timerResize("mat.resizeNonZeros()");
    mat.resizeNonZeros(totalNnz);
    timerResize.finish();

    timerAux.finish();

    timer::commonTimer timerFilling("filling");

    int* outer = mat.outerIndexPtr();
    int* ind = mat.innerIndexPtr();
    double* values = mat.valuePtr();

    outer[0] = 0;
    for (int i = 1; i < rows + 1; ++i) {
        outer[i] = outer[i - 1] + outerIndexes[i];
    }

#pragma omp parallel for
    for (int i = 0; i < std::ssize(globalRowInfosMerged); ++i) {
        RowInfo<12 * 12>& rowInfo = globalRowInfosMerged[i];
        const int row = rowInfo.row;
        const int start = outer[row];
        const int size = outer[row + 1] - start;

        for (int j = 0; j < size; ++j) {
            values[start + j] = rowInfo.values[j];
        }
        for (int j = 0; j < size; ++j) {
            ind[start + j] = rowInfo.columns[j];
        }
    }
}

void Mesh::stencil_Lmat2_OPT3(Operator& mat, const Domain& domain) const {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    // static int check_count = 0;
    // std::vector<Triplet> trips;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;
    // const int num_threads = std::min(omp_get_max_threads(), 128);

    const int rows = mat.rows();
    const int nthr = omp_get_max_threads();
    // const int nthr = 1;

    std::vector<std::vector<RowInfo<12>>> rowInfosTest(nthr);

    std::vector<uint8_t> nonZeroBlocks(rows);
#pragma omp parallel for
    for (int i = 0; i < rows; ++i) {
        nonZeroBlocks[i] = 0;
    }

    timer::commonTimer timerCountBlocks("count nnz blocks");
#pragma omp parallel num_threads(nthr)
    {
        timer::commonTimer timerOmp("OMP loop");
#pragma omp for collapse(3) schedule(dynamic, 8)
        for (int i = BORDER; i < max_i; ++i) {
            for (int j = BORDER; j < max_j; ++j) {
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = LmatX2[sind(i, j, k)];

                    // X component
                    computeNonZeroBlocksInRows<XIndexer, XIndexer, 0>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);
                    computeNonZeroBlocksInRows<XIndexer, YIndexer, 1>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);
                    computeNonZeroBlocksInRows<XIndexer, ZIndexer, 2>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);

                    // Y component
                    computeNonZeroBlocksInRows<YIndexer, XIndexer, 3>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);
                    computeNonZeroBlocksInRows<YIndexer, YIndexer, 4>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);
                    computeNonZeroBlocksInRows<YIndexer, ZIndexer, 5>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);

                    // Z component
                    computeNonZeroBlocksInRows<ZIndexer, XIndexer, 6>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);
                    computeNonZeroBlocksInRows<ZIndexer, YIndexer, 7>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);
                    computeNonZeroBlocksInRows<ZIndexer, ZIndexer, 8>(i, j, k, block, xSize, ySize, zSize, TOL,
                                                                      nonZeroBlocks);
                }
            }
        }
    }
    timerCountBlocks.finish();

    std::vector<int> nonZeroBlocksOuter(rows + 1);
    nonZeroBlocksOuter[0] = 0;
    int nonZeroBlocksCount = 0;
    for (int i = 0; i < rows; ++i) {
        nonZeroBlocksCount += nonZeroBlocks[i];
        nonZeroBlocks[i] = 0;
        nonZeroBlocksOuter[i + 1] = nonZeroBlocksCount;
    }
    std::cout << "total nnz blocks: " << nonZeroBlocksCount << std::endl;
    std::vector<RowInfo<12>> testBlocks(nonZeroBlocksCount);

    timer::commonTimer timerUnpackingNew("unpacking new");
#pragma omp parallel num_threads(nthr)
    {
        timer::commonTimer timerOmp("OMP loop");
#pragma omp for collapse(3) schedule(dynamic, 8)
        for (int i = BORDER; i < max_i; ++i) {
            for (int j = BORDER; j < max_j; ++j) {
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = LmatX2[sind(i, j, k)];

                    // X component
                    processComponent3<XIndexer, XIndexer, 0>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);
                    processComponent3<XIndexer, YIndexer, 1>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);
                    processComponent3<XIndexer, ZIndexer, 2>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);

                    // Y component
                    processComponent3<YIndexer, XIndexer, 3>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);
                    processComponent3<YIndexer, YIndexer, 4>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);
                    processComponent3<YIndexer, ZIndexer, 5>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);

                    // Z component
                    processComponent3<ZIndexer, XIndexer, 6>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);
                    processComponent3<ZIndexer, YIndexer, 7>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);
                    processComponent3<ZIndexer, ZIndexer, 8>(i, j, k, block, xSize, ySize, zSize, TOL, testBlocks,
                                                             nonZeroBlocks, nonZeroBlocksOuter);
                }
            }
        }
    }
    timerUnpackingNew.finish();

    timer::commonTimer timerMerge("merge");

    std::vector<int> blocksStarts;
    blocksStarts.reserve(std::ssize(testBlocks));
    blocksStarts.push_back(0);

    for (int i = 0, j = 0; i != std::ssize(testBlocks); i = j) {
        while (j < std::ssize(testBlocks) && testBlocks[i].row == testBlocks[j].row) {
            j += 1;
        }
        if (j != i) {
            blocksStarts.push_back(j);
        }
    }

    std::vector<RowInfo<12 * 12>> globalRowInfosMerged(blocksStarts.size() - 1);

#pragma omp parallel for
    for (int i = 0; i < std::ssize(blocksStarts) - 1; ++i) {
        globalRowInfosMerged[i].mergeFromOthers(blocksStarts[i + 1] - blocksStarts[i], &testBlocks[blocksStarts[i]]);
    }
    timerMerge.finish();

    timer::commonTimer timerAux("aux");
    int totalNnz = 0;
    std::vector<int> outerIndexes(rows + 1);
#pragma omp parallel for
    for (int i = 0; i < rows + 1; ++i) {
        outerIndexes[i] = 0;
    }

#pragma omp parallel for reduction(+ : totalNnz)
    for (int j = 0; j < std::ssize(globalRowInfosMerged); ++j) {
        outerIndexes[globalRowInfosMerged[j].row + 1] = globalRowInfosMerged[j].nnz;
        totalNnz += globalRowInfosMerged[j].nnz;
    }

    timer::commonTimer timerResize("mat.resizeNonZeros()");
    mat.resizeNonZeros(totalNnz);
    timerResize.finish();

    timerAux.finish();

    timer::commonTimer timerFilling("filling");

    int* outer = mat.outerIndexPtr();
    int* ind = mat.innerIndexPtr();
    double* values = mat.valuePtr();

    outer[0] = 0;
    for (int i = 1; i < rows + 1; ++i) {
        outer[i] = outer[i - 1] + outerIndexes[i];
    }

#pragma omp parallel for
    for (int i = 0; i < std::ssize(globalRowInfosMerged); ++i) {
        RowInfo<12 * 12>& rowInfo = globalRowInfosMerged[i];
        const int row = rowInfo.row;
        const int start = outer[row];
        const int size = outer[row + 1] - start;

        for (int j = 0; j < size; ++j) {
            values[start + j] = rowInfo.values[j];
        }
        for (int j = 0; j < size; ++j) {
            ind[start + j] = rowInfo.columns[j];
        }
    }
}

void Mesh::stencil_Lmat2_OPT4(Operator& mat, const Domain& domain, StencilLmat2_OPT4Workspace& workspace) const {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    // static int check_count = 0;
    // std::vector<Triplet> trips;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;
    // const int num_threads = std::min(omp_get_max_threads(), 128);

    const int rows = mat.rows();
    const int nthr = omp_get_max_threads();
    // const int nthr = 1;

    std::vector<std::vector<RowInfo<12>>> rowInfosTest(nthr);

    timer::commonTimer timerUnpackingNew("unpacking new");
#pragma omp parallel num_threads(nthr)
    {
        std::vector<RowInfo<12>>& localRowInfos = rowInfosTest[omp_get_thread_num()];
        localRowInfos.reserve(1024 * 1024 * 4);
        timer::commonTimer timerOmp("OMP loop");
#pragma omp for collapse(3) schedule(dynamic, 8)
        for (int i = BORDER; i < max_i; ++i) {
            for (int j = BORDER; j < max_j; ++j) {
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = LmatX2[sind(i, j, k)];

                    // X component
                    processComponent2<XIndexer, XIndexer, 0>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<XIndexer, YIndexer, 1>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<XIndexer, ZIndexer, 2>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);

                    // Y component
                    processComponent2<YIndexer, XIndexer, 3>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<YIndexer, YIndexer, 4>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<YIndexer, ZIndexer, 5>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);

                    // Z component
                    processComponent2<ZIndexer, XIndexer, 6>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<ZIndexer, YIndexer, 7>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                    processComponent2<ZIndexer, ZIndexer, 8>(i, j, k, block, xSize, ySize, zSize, TOL, localRowInfos);
                }
            }
        }
    }
    timerUnpackingNew.finish();

    timer::commonTimer timerPseudoSort("pseudo sort");

    int totalUnsortedRows = 0;
    for (int i = 0; i < nthr; ++i) {
        totalUnsortedRows += std::ssize(rowInfosTest[i]);
    }

    timer::commonTimer timerInitForSort("init vectors of ints");
    // Indexes of each arrays is: row, orig thread owner and orig. number in its buffer
    std::vector<std::array<int, 3>>& rowPositionsUnsorted = workspace.rowPositionsUnsorted;
    std::vector<std::array<int, 3>>& rowPositionsSorted = workspace.rowPositionsSorted;
    if (rowPositionsSorted.capacity() < static_cast<size_t>(totalUnsortedRows)) {
        // Avoid memory reallocation, since old data is not necessary here
        std::vector<std::array<int, 3>> tmp1;
        std::vector<std::array<int, 3>> tmp2;
        tmp1.reserve(totalUnsortedRows * 2);
        tmp2.reserve(totalUnsortedRows * 2);
        rowPositionsUnsorted = std::move(tmp1);
        rowPositionsSorted = std::move(tmp2);
    }
    rowPositionsUnsorted.resize(totalUnsortedRows);
    rowPositionsSorted.resize(totalUnsortedRows);
    timerInitForSort.finish();

    timer::commonTimer timerSetUnsortedRowPos("set row positions unsorted");

    std::vector<int> offsets(nthr);
    offsets[0] = 0;
    for (int i = 1; i < nthr; ++i) {
        offsets[i] = std::ssize(rowInfosTest[i - 1]) + offsets[i - 1];
    }
#pragma omp parallel for num_threads(nthr)
    for (int i = 0; i < nthr; ++i) {
        const std::vector<RowInfo<12>>& localRowInfos = rowInfosTest[i];

        for (int j = 0; j < std::ssize(localRowInfos); ++j) {
            rowPositionsUnsorted[offsets[i] + j] = {localRowInfos[j].row, i, j};
        }
    }

    timerSetUnsortedRowPos.finish();

    timer::commonTimer timerComputingNnzBlocks("compute nnz block count in rows");

    std::vector<uint8_t> nonZeroBlocks(rows);
#pragma omp parallel for
    for (int i = 0; i < rows; ++i) {
        nonZeroBlocks[i] = 0;
    }

#pragma omp parallel for
    for (int i = 0; i < nthr; ++i) {
        const std::vector<RowInfo<12>>& localRowInfos = rowInfosTest[i];
        for (int j = 0; j < std::ssize(localRowInfos); ++j) {
            std::atomic_ref<uint8_t> toUpd(nonZeroBlocks[localRowInfos[j].row]);
            toUpd += 1;
        }
    }

    timerComputingNnzBlocks.finish();

    timer::commonTimer timerMulti("set blocks bound and zero aux array");

    std::vector<int> nonZeroBlocksOuter(rows + 1);
    nonZeroBlocksOuter[0] = 0;
    int nonZeroBlocksCount = 0;
    for (int i = 0; i < rows; ++i) {
        nonZeroBlocksCount += nonZeroBlocks[i];
        nonZeroBlocks[i] = 0;
        nonZeroBlocksOuter[i + 1] = nonZeroBlocksCount;
    }

    timerMulti.finish();

    timer::commonTimer timerSetRowPosSorted("set row positions sorted");
#pragma omp parallel for
    for (int i = 0; i < totalUnsortedRows; ++i) {
        const int row = rowPositionsUnsorted[i][0];
        std::atomic_ref<uint8_t> blockIxRef(nonZeroBlocks[row]);
        const int blockIx = (blockIxRef++) + nonZeroBlocksOuter[row];
        rowPositionsSorted[blockIx] = rowPositionsUnsorted[i];
    }

    timerSetRowPosSorted.finish();

    timer::commonTimer timerSetBlocksBounds("set block bounds in glob. array");

    std::vector<int> blocksStarts;
    blocksStarts.reserve(std::ssize(rowPositionsSorted));
    blocksStarts.push_back(0);

    for (int i = 0, j = 0; i != std::ssize(rowPositionsSorted); i = j) {
        while (j < std::ssize(rowPositionsSorted) && rowPositionsSorted[i][0] == rowPositionsSorted[j][0]) {
            j += 1;
        }
        if (j != i) {
            blocksStarts.push_back(j);
        }
    }

    timerSetBlocksBounds.finish();

    timerPseudoSort.finish();

    timer::commonTimer timerMerge("merge");

    std::vector<RowInfo<12 * 12>> globalRowInfosMerged(blocksStarts.size() - 1);

#pragma omp parallel
    {
        timer::commonTimer timerOmp("OMP section");
        constexpr int maxMergedBlocks = 12 * 12;
        std::array<RowInfo<12>, maxMergedBlocks> tmpStorage;
#pragma omp for
        for (int i = 0; i < std::ssize(blocksStarts) - 1; ++i) {
            const int start = blocksStarts[i];
            const int end = blocksStarts[i + 1];
            assert(end - start <= maxMergedBlocks);
            for (int j = start; j < end; ++j) {
                const auto [row, thread, index] = rowPositionsSorted[j];
                tmpStorage[j - start] = rowInfosTest[thread][index];
            }
            globalRowInfosMerged[i].mergeFromOthers(end - start, &tmpStorage[0]);
        }
    }
    timerMerge.finish();

    timer::commonTimer timerAux("aux");
    int totalNnz = 0;
    std::vector<int> outerIndexes(rows + 1);
#pragma omp parallel for
    for (int i = 0; i < rows + 1; ++i) {
        outerIndexes[i] = 0;
    }

#pragma omp parallel for reduction(+ : totalNnz)
    for (int i = 0; i < std::ssize(globalRowInfosMerged); ++i) {
        outerIndexes[globalRowInfosMerged[i].row + 1] = globalRowInfosMerged[i].nnz;
        totalNnz += globalRowInfosMerged[i].nnz;
    }

    timer::commonTimer timerResize("mat.resizeNonZeros()", totalNnz);
    /* Usage of timerResize.finish();
     * leads to memory reallocation, i.e. it calls realloc inside Eigen code. In order to avoid
     * memory single-thread copy inside Eigen at least in some cases (and do not modify it's source code), reserve
     * memory for non-zero indexes with factor 2 here.
     */
    if (mat.data().allocatedSize() < totalNnz) {
        timer::commonTimer timerReserve("realloc in eigen buffer");
        // reserve() acquires increases buffer BY value, not a TO value
        mat.data().reserve(totalNnz * 2 - mat.data().allocatedSize());
    }
    mat.resizeNonZeros(totalNnz);
    timerResize.finish();

    timerAux.finish();

    timer::commonTimer timerFilling("filling");

    int* outer = mat.outerIndexPtr();
    int* ind = mat.innerIndexPtr();
    double* values = mat.valuePtr();

    outer[0] = 0;
    for (int i = 1; i < rows + 1; ++i) {
        outer[i] = outer[i - 1] + outerIndexes[i];
    }

#pragma omp parallel for
    for (int i = 0; i < std::ssize(globalRowInfosMerged); ++i) {
        RowInfo<12 * 12>& rowInfo = globalRowInfosMerged[i];
        const int row = rowInfo.row;
        const int start = outer[row];
        const int size = outer[row + 1] - start;

        for (int j = 0; j < size; ++j) {
            values[start + j] = rowInfo.values[j];
        }
        for (int j = 0; j < size; ++j) {
            ind[start + j] = rowInfo.columns[j];
        }
    }
}

void Mesh::stencil_Lmat2(Operator& mat, const Domain& domain) const {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    static int check_count = 0;
    // std::vector<Triplet> trips;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;
    const int num_threads = std::min(omp_get_max_threads(), 128);

    // Оценка количества итераций и элементов для каждого потока.
    // Если общее число итераций равно:
    //   (max_i - BORDER) * (max_j - BORDER) * (max_k - BORDER)
    // и для каждой итерации добавляется 9 элементов (3x3 окрестность),
    // то можно примерно вычислить estimated_per_thread.
    size_t total_iterations = static_cast<size_t>(max_i - BORDER) * (max_j - BORDER) * (max_k - BORDER);
    size_t estimated_per_thread = (total_iterations / num_threads) * 9;

    // Локальные векторы для каждого потока с предварительным резервированием
    // памяти.
    std::vector<std::vector<Triplet>> local_vectors(num_threads);

    timer::commonTimer timer1("parallel sections");

#pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        // First touch - выделяем память в том же потоке, который будет её
        // использовать
        local_vectors[tid].reserve(estimated_per_thread);
    }
#pragma omp parallel num_threads(num_threads)
    {
        timer::flatTimer timerLocal("main OMP section");

        int tid = omp_get_thread_num();

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const Block& block = LmatX2[sind(i, j, k)];

                    // if (i != 230 || j != 202 || k != 1) {
                    //     continue;
                    // }

                    // std::cout << "################################### ref proceeding block " << i << " " << j << " "
                    //           << k << std::endl;

                    // X component
                    processComponent<XIndexer, XIndexer, 0>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<XIndexer, YIndexer, 1>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<XIndexer, ZIndexer, 2>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);

                    // Y component
                    processComponent<YIndexer, XIndexer, 3>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<YIndexer, YIndexer, 4>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<YIndexer, ZIndexer, 5>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);

                    // Z component
                    processComponent<ZIndexer, XIndexer, 6>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<ZIndexer, YIndexer, 7>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                    processComponent<ZIndexer, ZIndexer, 8>(i, j, k, block, local_vectors[tid], xSize, ySize, zSize,
                                                            TOL);
                }

        // Сортировка локального вектора по (row, col)
        auto& vec = local_vectors[tid];
        timer::flatTimer timerDuplications("remove duplications", vec.size());
        timer::flatTimer timerSort("sort vec", vec.size());
        std::sort(vec.begin(), vec.end(), compareTriplets);
        timerSort.finish();

        // Устранение дубликатов: проход по отсортированному вектору, складываем
        // значения для одинаковых ключей
        if (!vec.empty()) {
            size_t index = 0;
            for (size_t j = 1; j < vec.size(); ++j) {
                if (vec[index].row() == vec[j].row() && vec[index].col() == vec[j].col()) {
                    vec[index].value() += vec[j].value();
                } else {
                    ++index;
                    vec[index] = vec[j];
                }
            }
            vec.resize(index + 1);
        }
    }

    timer1.finish();
    timer::commonTimer timer2("section 2");

    pmms::PMMSOptions opt;
    opt.useSampling = true;
    opt.usePWayMerge = true;
    opt.oversample = 4;
    std::vector<Triplet> test = pmms::parallelMultiwayMergeSort<Triplet>(std::move(local_vectors), num_threads, opt);
    if (!test.empty()) {
        size_t index = 0;
        for (size_t j = 1; j < test.size(); ++j) {
            if (test[index].row() == test[j].row() && test[index].col() == test[j].col()) {
                test[index].value() += test[j].value();
            } else {
                ++index;
                test[index] = test[j];
            }
        }
        test.resize(index + 1);
    }

    timer2.finish();
    timer::commonTimer timer3("section 3");

    check_count++;
    // В non_empty[0] теперь находится глобальный вектор, уже
    // отсортированный и с устранёнными дубликатами.
    if (check_count < 0) {
        const std::vector<Triplet> trips = multyPhaseMerge(local_vectors);
        if (!equalVecsTriplets(test, trips))
            LOG_STEP("Not Equal!\n");
    }

    // std::vector<Trip> trips2;
    // trips2.reserve(trips.size());
    // for (long i = 0; i < trips.size(); i++){
    //     trips2.emplace_back(trips[i].row, trips[i].col, trips[i].value);
    // }

    // mat.setFromTriplets(trips.begin(), trips.end());

    timer3.finish();

    optimizedSetFromSortedTriplets(mat, test);

    LOG_STEP("Matrix L (block) was created." << " trips size: " << test.size() << std::endl);
}

// TODO implement parallelMultiwayMergeSort and merge to convert block to csr
// (it is dublicated)
void Mesh::stencil_Lmat2_NGP(Operator& mat, const Domain& domain) {
    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;

    // std::vector<Triplet> trips;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;
    const int num_threads = std::min(omp_get_max_threads(), 128);

    // Оценка количества итераций и элементов для каждого потока.
    // Если общее число итераций равно:
    //   (max_i - BORDER) * (max_j - BORDER) * (max_k - BORDER)
    // и для каждой итерации добавляется 9 элементов (3x3 окрестность),
    // то можно примерно вычислить estimated_per_thread.
    size_t total_iterations = static_cast<size_t>(max_i - BORDER) * (max_j - BORDER) * (max_k - BORDER);
    size_t estimated_per_thread = (total_iterations / num_threads) * 9;

    // Локальные векторы для каждого потока с предварительным резервированием
    // памяти.
    std::vector<std::vector<Triplet>> local_vectors(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        local_vectors[t].reserve(estimated_per_thread);
    }

#pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();

        std::vector<Triplet> local_trips;
        local_trips.reserve(16 * 9 * (max_i - BORDER) / omp_get_num_threads());

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX_NGP.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = LmatX_NGP[sind(i, j, k)];

                    // X component
                    processComponent<XIndexerNGP, XIndexerNGP, 0>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<XIndexerNGP, YIndexerNGP, 1>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<XIndexerNGP, ZIndexerNGP, 2>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);

                    // Y component
                    processComponent<YIndexerNGP, XIndexerNGP, 3>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<YIndexerNGP, YIndexerNGP, 4>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<YIndexerNGP, ZIndexerNGP, 5>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);

                    // Z component
                    processComponent<ZIndexerNGP, XIndexerNGP, 6>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<ZIndexerNGP, YIndexerNGP, 7>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                    processComponent<ZIndexerNGP, ZIndexerNGP, 8>(i, j, k, block, local_vectors[tid], xSize, ySize,
                                                                  zSize, TOL);
                }

        // Сортировка локального вектора по (row, col)
        auto& vec = local_vectors[tid];
        std::sort(vec.begin(), vec.end(), compareTriplets);

        // Устранение дубликатов: проход по отсортированному вектору, складываем
        // значения для одинаковых ключей
        if (!vec.empty()) {
            size_t index = 0;
            for (size_t j = 1; j < vec.size(); ++j) {
                if (vec[index].row() == vec[j].row() && vec[index].col() == vec[j].col()) {
                    vec[index].value() += vec[j].value();
                } else {
                    ++index;
                    vec[index] = vec[j];
                }
            }
            vec.resize(index + 1);
        }
    }

    // Собираем только непустые локальные векторы для дальнейшего слияния
    std::vector<std::vector<Triplet>> non_empty;
    for (auto& v : local_vectors) {
        if (!v.empty()) {
            non_empty.push_back(std::move(v));
        }
    }

    // Многофазное слияние: объединяем пары отсортированных векторов параллельно
    // с предварительным резервированием памяти
    while (non_empty.size() > 1) {
        size_t new_size = (non_empty.size() + 1) / 2;
        std::vector<std::vector<Triplet>> new_vectors(new_size);
        size_t pairs = non_empty.size() / 2;

#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < pairs; ++i) {
            const auto& left = non_empty[2 * i];
            const auto& right = non_empty[2 * i + 1];
            // Резервируем память для слияния двух векторов
            size_t merged_capacity = left.size() + right.size();
            std::vector<Triplet> merged;
            merged.reserve(merged_capacity);

            size_t li = 0, ri = 0;
            while (li < left.size() && ri < right.size()) {
                if (compareTriplets(left[li], right[ri])) {
                    merged.push_back(left[li]);
                    ++li;
                } else if (compareTriplets(right[ri], left[li])) {
                    merged.push_back(right[ri]);
                    ++ri;
                } else {
                    // Если ключи равны – складываем значения
                    Triplet t = left[li];
                    t.value() += right[ri].value();
                    merged.push_back(t);
                    ++li;
                    ++ri;
                }
            }
            while (li < left.size()) {
                merged.push_back(left[li++]);
            }
            while (ri < right.size()) {
                merged.push_back(right[ri++]);
            }
            new_vectors[i] = std::move(merged);
        }
        // Если число векторов нечётное – последний переносим без изменений
        if (non_empty.size() % 2 == 1) {
            new_vectors.back() = std::move(non_empty.back());
        }
        non_empty = std::move(new_vectors);
    }

    // В non_empty[0] теперь находится глобальный вектор, уже отсортированный и
    // с устранёнными дубликатами.
    const std::vector<Triplet>& trips = non_empty.empty() ? std::vector<Triplet>() : non_empty.front();

    // std::vector<Trip> trips2;
    // trips2.reserve(trips.size());
    // for (long i = 0; i < trips.size(); i++){
    //     trips2.emplace_back(trips[i].row, trips[i].col, trips[i].value);
    // }

    // mat.setFromTriplets(trips.begin(), trips.end());
    optimizedSetFromSortedTriplets(mat, trips);
    LOG_STEP("Matrix L (block) was created." << " trips size: " << trips.size() << std::endl);
}

void Mesh::stencil_curlB(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    auto addRowIfInside = [&](int i, int j, int k, int comp) -> bool {
        return domain.is_inside_node_periodic(i, j, k, FieldType::ELECTRIC, comp, is_periodic);
    };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int im = bc_handler.wrap_index(i - 1, 0, domain);
                const int jm = bc_handler.wrap_index(j - 1, 1, domain);
                const int km = bc_handler.wrap_index(k - 1, 2, domain);

                const int vindx = domain.vind(i, j, k, 0);   // Ex
                const int vindy = domain.vind(i, j, k, 1);   // Ey
                const int vindz = domain.vind(i, j, k, 2);   // Ez

                // (x)[i+1/2,j,k]
                if (addRowIfInside(i, j, k, 0)) {
                    // ( Bz[i+1/2,j+1/2,k] - Bz[i+1/2,j-1/2,k] ) / dy
                    double val = 1.0 / dy;
                    trips.emplace_back(vindx, domain.vind(i, j, k, 2), val);
                    trips.emplace_back(vindx, domain.vind(i, jm, k, 2), -val);
                    // - ( By[i+1/2,j,k+1/2] - By[i+1/2,j,k-1/2] ) / dz
                    val = -1.0 / dz;
                    trips.emplace_back(vindx, domain.vind(i, j, k, 1), val);
                    trips.emplace_back(vindx, domain.vind(i, j, km, 1), -val);
                }
                // (y)[i,j+1/2,k]
                if (addRowIfInside(i, j, k, 1)) {
                    // ( Bx[i,j+1/2,k+1/2] - Bx[i,j+1/2,k-1/2] ) / dz
                    double val = 1.0 / dz;
                    trips.emplace_back(vindy, domain.vind(i, j, k, 0), val);
                    trips.emplace_back(vindy, domain.vind(i, j, km, 0), -val);
                    // -( Bz[i+1/2,j+1/2,k] - Bz[i-1/2,j+1/2,k] ) / dx
                    val = -1.0 / dx;
                    trips.emplace_back(vindy, domain.vind(i, j, k, 2), val);
                    trips.emplace_back(vindy, domain.vind(im, j, k, 2), -val);
                }
                // (z)[i,j,k+1/2]
                if (addRowIfInside(i, j, k, 2)) {
                    // ( By[i+1/2,j,k+1/2] - By[i-1/2,j,k+1/2] ) / dx
                    double val = 1.0 / dx;
                    trips.emplace_back(vindz, domain.vind(i, j, k, 1), val);
                    trips.emplace_back(vindz, domain.vind(im, j, k, 1), -val);
                    // -( Bx[i,j+1/2,k+1/2] - Bx[i,j-1/2,k+1/2] ) / dy
                    val = -1.0 / dy;
                    trips.emplace_back(vindz, domain.vind(i, j, k, 0), val);
                    trips.emplace_back(vindz, domain.vind(i, jm, k, 0), -val);
                }
                //  bc_handler.modify_curlB_stencil(i, j, k, trips, domain);
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_curlE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 12;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    // Проверяет, принадлежит ли магнитный узел области (с учётом периодичности)
    auto addRowIfInside = [&](int i, int j, int k, int comp) -> bool {
        return domain.is_inside_node_periodic(i, j, k, FieldType::MAGNETIC, comp, is_periodic);
    };

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                const int ip = bc_handler.wrap_index(i + 1, 0, domain);
                const int jp = bc_handler.wrap_index(j + 1, 1, domain);
                const int kp = bc_handler.wrap_index(k + 1, 2, domain);

                const int vindx = domain.vind(i, j, k, 0);   // Bx
                const int vindy = domain.vind(i, j, k, 1);   // By
                const int vindz = domain.vind(i, j, k, 2);   // Bz

                // (x)[i,j+1/2,k+1/2]
                if (addRowIfInside(i, j, k, 0)) {
                    // ( Ez[i,j+1,k+1/2] - Ez[i,j,k+1/2] ) / dy
                    double val = 1.0 / dy;
                    trips.emplace_back(vindx, domain.vind(i, jp, k, 2), val);
                    trips.emplace_back(vindx, domain.vind(i, j, k, 2), -val);
                    // - ( Ey[i,j+1/2,k+1] - Ey[i,j+1/2,k] ) / dz
                    val = -1.0 / dz;
                    trips.emplace_back(vindx, domain.vind(i, j, kp, 1), val);
                    trips.emplace_back(vindx, domain.vind(i, j, k, 1), -val);
                }

                // (y)[i+1/2,j,k+1/2]
                if (addRowIfInside(i, j, k, 1)) {
                    // ( Ex[i+1/2,j,k+1] - Ex[i+1/2,j,k] ) / dz
                    double val = 1.0 / dz;
                    trips.emplace_back(vindy, domain.vind(i, j, kp, 0), val);
                    trips.emplace_back(vindy, domain.vind(i, j, k, 0), -val);
                    // - ( Ez[i+1,j,k+1/2] - Ez[i,j,k+1/2] ) / dx
                    val = -1.0 / dx;
                    trips.emplace_back(vindy, domain.vind(ip, j, k, 2), val);
                    trips.emplace_back(vindy, domain.vind(i, j, k, 2), -val);
                }

                // (z)[i+1/2,j+1/2,k]
                if (addRowIfInside(i, j, k, 2)) {
                    // ( Ey[i+1,j+1/2,k] - Ey[i,j+1/2,k] ) / dx
                    double val = 1.0 / dx;
                    trips.emplace_back(vindz, domain.vind(ip, j, k, 1), val);
                    trips.emplace_back(vindz, domain.vind(i, j, k, 1), -val);
                    // - ( Ex[i+1/2,j+1,k] - Ex[i+1/2,j,k] ) / dy
                    val = -1.0 / dy;
                    trips.emplace_back(vindz, domain.vind(i, jp, k, 0), val);
                    trips.emplace_back(vindz, domain.vind(i, j, k, 0), -val);
                }

                //  bc_handler.modify_curlE_stencil(i, j, k, trips, domain);
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}

void Mesh::stencil_divE(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    std::vector<Trip> trips;
    const auto size = domain.size();
    int totalSize = size.x() * size.y() * size.z() * 6;
    trips.reserve(totalSize);

    const double dx = domain.cell_size().x();
    const double dy = domain.cell_size().y();
    const double dz = domain.cell_size().z();

    bool is_periodic[3] = {bc_handler.is_periodic(0), bc_handler.is_periodic(1), bc_handler.is_periodic(2)};

    for (int i = 0; i < size.x(); i++) {
        for (int j = 0; j < size.y(); j++) {
            for (int k = 0; k < size.z(); k++) {
                if (!domain.is_inside_node_periodic(i, j, k, FieldType::DENSITY, 0, is_periodic))
                    continue;

                const int im = bc_handler.wrap_index(i - 1, 0, domain);
                const int jm = bc_handler.wrap_index(j - 1, 1, domain);
                const int km = bc_handler.wrap_index(k - 1, 2, domain);

                const int sindx = domain.sind(i, j, k);

                // [i,j,k]
                // ( Ex[i+1/2,j,k] - Ex[i-1/2,j,k] ) / dx
                double val = 1.0 / dx;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 0), val));
                trips.push_back(Trip(sindx, domain.vind(im, j, k, 0), -val));
                // ( Ey[i,j+1/2,k] - Ey[i,j-1/2,k] ) / dy
                val = 1.0 / dy;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 1), val));
                trips.push_back(Trip(sindx, domain.vind(i, jm, k, 1), -val));
                // ( Ez[i,j,k+1/2] - Ez[i,j,k-1/2] ) / dz
                val = 1.0 / dz;
                trips.push_back(Trip(sindx, domain.vind(i, j, k, 2), val));
                trips.push_back(Trip(sindx, domain.vind(i, j, km, 2), -val));
                // TODO: add bc_handler modify stencil
            }
        }
    }
    mat.setFromTriplets(trips.begin(), trips.end());
}
void Mesh::stencil_smooth_1d(Operator& mat, const Domain& domain, int dim) {
    constexpr int COMPONENTS = 3;
    constexpr int STENCIL = 3;
    const auto size = domain.size();
    const int nx = size.x();
    const int ny = size.y();
    const int nz = size.z();

    const int totalCells = nx * ny * nz;
    // checker
    mat.resize(domain.total_size() * 3, domain.total_size() * 3);

    if (dim != Axis::X && dim != Axis::Y && dim != Axis::Z) {
        std::cout << "Error: Invalid dimension." << std::endl;
        return;
    }
    std::vector<Trip> trips;
    trips.reserve(static_cast<size_t>(totalCells * COMPONENTS * STENCIL));

    auto addTriplet = [](std::vector<Trip>& trips, const Domain& domain, int row, int col, double val) {
        bool onArea = domain.is_inside_node(row, FieldType::ELECTRIC);
        if (onArea) {
            trips.push_back(Trip(row, col, val));
        }
    };

    // auto imod = [](int a, int m) {
    //     // % m - 1
    //     if (a > m - 1) return 3;
    //     if (a < 0) return m - 4;
    //     return a;
    // };

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                int im = i - 1, ip = i + 1;
                int jm = j - 1, jp = j + 1;
                int km = k - 1, kp = k + 1;

                // if (dim == Axis::X && domain.is_periodic_bound(X)) {
                //     im = imod(i - 1, nx);
                //     ip = imod(i + 1, nx);
                // } else if (dim == Axis::Y && domain.is_periodic_bound(Y)) {
                //     jm = imod(j - 1, ny);
                //     jp = imod(j + 1, ny);
                // } else if (dim == Axis::Z &&
                //            domain.is_periodic_bound(Z)) {
                //     km = imod(k - 1, nz);
                //     kp = imod(k + 1, nz);
                // }

                for (int ax = 0; ax < COMPONENTS; ++ax) {
                    const int row = vind(i, j, k, ax);

                    int col_m, col_c, col_p;
                    if (dim == Axis::X) {
                        col_m = vind(im, j, k, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(ip, j, k, ax);
                    } else if (dim == Axis::Y) {
                        col_m = vind(i, jm, k, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(i, jp, k, ax);
                    } else {
                        col_m = vind(i, j, km, ax);
                        col_c = vind(i, j, k, ax);
                        col_p = vind(i, j, kp, ax);
                    }

                    addTriplet(trips, domain, row, col_m, 0.25);
                    addTriplet(trips, domain, row, col_c, 0.5);
                    addTriplet(trips, domain, row, col_p, 0.25);
                }
            }
        }
    }

    mat.setFromTriplets(trips.begin(), trips.end());
}
