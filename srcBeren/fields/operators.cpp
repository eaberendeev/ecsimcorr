#include "operators.h"

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
    RowInfo(int rowIn) : row(rowIn) {
    }
    RowInfo() {
    }

    std::array<double, maxNnz> values;
    std::array<int, maxNnz> columns;
    int row;
    int nnz{0};

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
        int its[count]{0};
        int counter = 0;

        for (int i = 1; i < count; ++i) {
            assert(others[i].row == others[0].row);
        }

        row = others[0].row;

        while (true) {
            int smallestCol = std::numeric_limits<int>::max();
            for (int i = 0; i < count; ++i) {
                if (its[i] < others[i].nnz)
                    smallestCol = std::min(smallestCol, others[i].columns[its[i]]);
            }

            if (smallestCol == std::numeric_limits<int>::max()) {
                break;
            }

            double acc = 0.0;
            for (int i = 0; i < count; ++i) {
                if (its[i] < others[i].nnz && others[i].columns[its[i]] == smallestCol) {
                    acc += others[i].values[its[i]];
                    its[i] += 1;
                }
            }

            assert(counter < maxNnz);
            columns[counter] = smallestCol;
            values[counter] = acc;
            counter += 1;
            nnz = counter;
        }
    }
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
}

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
void processComponent2(int i_cell, int j_cell, int k_cell, const Block_t& block, [[maybe_unused]] int Nx, int Ny,
                       int Nz, double tolerance, std::vector<RowInfo<12 * 9>>& rowInfos) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1) {
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1) {
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row = vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                     k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                const int rowIdx = RowIdx::calculate(x1, y1, z1);
                assert(rowIdx >= 0 && rowIdx < 12);
                RowInfo<12 * 9> rowInfo(row);

                for (int x2 = 0; x2 < ColIdx::size_x; ++x2) {
                    for (int y2 = 0; y2 < ColIdx::size_y; ++y2) {
                        for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                            const double val = block(rowIdx, ColIdx::calculate(x2, y2, z2), DIR);

                            // std::cout << "test coords and val: " << x1 << " " << y1 << " " << z1 << " " << x2 << " "
                            //           << y2 << " " << z2 << " " << val << " for row " << row << std::endl;

                            if (std::abs(val) > tolerance) {
                                // std::cout << "test add to row " << row << std::endl;
                                const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                                     k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                                rowInfo.addValue(col, val);
                            }
                        }
                    }
                }
                if (rowInfo.nnz != 0) {
                    // std::cout << "non zero element for row " << row << std::endl;
                    rowInfo.sort();
                    rowInfos.push_back(rowInfo);
                }
            }
        }
    }
}

void Mesh::stencil_Lmat2_OPT(Operator& mat, const Domain& domain) const {
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

    std::vector<std::vector<RowInfo<12 * 9>>> rowInfosTest(nthr);
    for (int i = 0; i < nthr; ++i) {
        rowInfosTest[i].reserve(1024);
    }

    timer::commonTimer timerUnpackingNew("unpacking new");
    // #pragma omp for collapse(3) schedule(dynamic, 8) nowait
    for (int i = BORDER; i < max_i; ++i) {
        for (int j = BORDER; j < max_j; ++j) {
            for (int k = BORDER; k < max_k; ++k) {
                if (!LmatX2.non_zeros[sind(i, j, k)])
                    continue;
                const auto& block = LmatX2[sind(i, j, k)];

                // if (i != 230 || j != 202 || k != 1) {
                //     continue;
                // }

                // std::cout << "################################### test proceeding block " << i << " " << j << " " <<
                // k
                //           << std::endl;

                std::vector<RowInfo<12 * 9>>& localRowInfos = rowInfosTest[omp_get_thread_num()];

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
    timerUnpackingNew.finish();

    timer::commonTimer timerSortNew("sort new");
    for (int i = 0; i < nthr; ++i) {
        std::vector<RowInfo<12 * 9>>& localRowInfos = rowInfosTest[i];

        std::sort(localRowInfos.begin(), localRowInfos.end(),
                  [](const RowInfo<12 * 9>& a, const RowInfo<12 * 9>& b) { return a.row < b.row; });
    }
    timerSortNew.finish();

    std::vector<std::vector<RowInfo<12 * 12 * 9>>> localRowInfosMerged(nthr);
    for (int i = 0; i < nthr; ++i) {
        localRowInfosMerged[i].reserve(1024);
    }
    timer::commonTimer timerMergeNew("merge New");
    for (int i = 0; i < nthr; ++i) {
        std::vector<RowInfo<12 * 9>>& localRowInfos = rowInfosTest[i];

        for (int i = 0, j = 0; i != std::ssize(localRowInfos); i = j) {
            while (localRowInfos[i].row == localRowInfos[j].row) {
                j += 1;
            }

            assert(j - i <= 12 * 9 && j - i > 0);

            RowInfo<12 * 12 * 9> rowInfo;
            rowInfo.mergeFromOthers(j - i, &localRowInfos[i]);

            localRowInfosMerged[0].push_back(rowInfo);
            // std::cout << "push line # " << localRowInfosMerged[0].size() << " from range " << i << " ... " << j
            //           << " for row " << rowInfo.row << std::endl;
        }
    }
    timerMergeNew.finish();

    // auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * ySize * zSize + j * zSize + k); };

    // timer::commonTimer timerCountingRows("getting non zero rows");
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

    // timer::commonTimer timerUnpacking("unpacking");
    // // #pragma omp parallel for num_threads(nthr)
    // for (int row = 0; row < rows; ++row) {
    //     // for (int i = 0; i < std::ssize(nonEmptyRows); ++i) {
    //     // const int row = nonEmptyRows[i];
    //     RowInfo rowInfo(row);

    //     // X component
    //     if (row % 3 == 0)
    //         processRow2<XIndexer, 0>(*this, row, max_i, max_j, max_k, BORDER, xSize, ySize, zSize, TOL, rowInfo);
    //     // Y component
    //     if (row % 3 == 1)
    //         processRow2<YIndexer, 3>(*this, row, max_i, max_j, max_k, BORDER, xSize, ySize, zSize, TOL, rowInfo);
    //     // Z component
    //     if (row % 3 == 2)
    //         processRow2<ZIndexer, 6>(*this, row, max_i, max_j, max_k, BORDER, xSize, ySize, zSize, TOL, rowInfo);

    //     if (rowInfo.nnz != 0) {
    //         rowInfo.sort();
    //         localRowInfos[omp_get_thread_num()].push_back(rowInfo);
    //     }
    // }
    // timerUnpacking.finish();

    int totalNnz = 0;
    std::vector<int> outerIndexes(rows + 1);
    std::fill(outerIndexes.begin(), outerIndexes.end(), 0);

    // #pragma omp parallel for reduction(+ : totalNnz)
    for (int i = 0; i < nthr; ++i) {
        const std::vector<RowInfo<>>& rowInfos = localRowInfosMerged[i];
        for (int j = 0; j < std::ssize(rowInfos); ++j) {
            outerIndexes[rowInfos[j].row + 1] = rowInfos[j].nnz;
            totalNnz += rowInfos[j].nnz;
        }
    }
    std::cout << "total Nnz: " << totalNnz << std::endl;

    for (int i = 1; i < rows + 1; ++i) {
        outerIndexes[i] += outerIndexes[i - 1];
    }

    mat.resizeNonZeros(totalNnz);
    int* outer = mat.outerIndexPtr();
    int* ind = mat.innerIndexPtr();
    // int* innerNnzRes = res.innerNonZeroPtr();

    double* values = mat.valuePtr();

    timer::commonTimer timerFilling("filling");

    outer[0] = 0;
#pragma omp parallel for
    for (int i = 1; i < rows + 1; ++i) {
        outer[i] = outerIndexes[i];
    }

    for (int i = 0; i < std::ssize(localRowInfosMerged[0]); ++i) {
        RowInfo<12 * 12 * 9>& rowInfo = localRowInfosMerged[0][i];
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
    if (check_count < 20) {
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
