#include <atomic>
#include <memory>
#include <source_location>

#include "Mesh.h"
#include "Shape.h"
#include "World.h"
#include "config.h"
#include "env_options.h"
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

template <int maxNnz = 12 * 12 * 9>
struct RowBlock {
    RowBlock(int rowIn) : row(rowIn), nnz(0) {
    }
    RowBlock() {
    }

    ~RowBlock() {
    }

    template <int otherNnz>
    RowBlock(int count, const RowBlock<otherNnz>* others) {
        mergeFromOthers(count, others);
    }

    RowBlock& operator=(const RowBlock<maxNnz>& other) {
        nnz = other.nnz;
        row = other.row;
        assert(nnz <= maxNnz);

        std::copy_n(other.values.begin(), nnz, values.begin());
        std::copy_n(other.columns.begin(), nnz, columns.begin());

        return *this;
    }

    RowBlock(const RowBlock& other) {
        nnz = other.nnz;
        row = other.row;

        std::copy_n(other.values.begin(), nnz, values.begin());
        std::copy_n(other.columns.begin(), nnz, columns.begin());
    }

    bool operator!=(const RowBlock& other) {
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

    void push_back_value(int col, double val) {
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
    void mergeFromOthers(int count, const RowBlock<otherNnz>* others) {
        int smallestCols[(maxNnz + otherNnz - 1) / otherNnz];

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
            int smallestColsCount = 0;
            double acc = 0;
            for (int i = 0; i < count; ++i) {
                if (its[i] < othersNnz[i]) {
                    const int otherCol = others[i].columns[its[i]];
                    if (otherCol == smallestCol) {
                        acc += others[i].values[its[i]];
                        smallestCols[smallestColsCount] = i;
                        smallestColsCount += 1;
                    } else if (otherCol < smallestCol) {
                        smallestCol = otherCol;
                        acc = others[i].values[its[i]];
                        smallestCols[0] = i;
                        smallestColsCount = 1;
                    }
                }
            }

            if (smallestCol == std::numeric_limits<int>::max()) {
                break;
            }

            for (int i = 0; i < smallestColsCount; ++i) {
                its[smallestCols[i]] += 1;
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

template <typename T>
struct SimpleArrayBuffer : public std::vector<T> {
    using std::vector<T>::size;
    using std::vector<T>::capacity;
    using std::vector<T>::reserve;

    T& operator()(const int64_t i) {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    const T& operator()(const int64_t i) const {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    T& operator[](const int64_t i) {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    const T& operator[](const int64_t i) const {
        assert(i >= 0 && i < std::ssize(*this));
        return std::vector<T>::operator[](i);
    }

    // resize buffer to fit new size, does not carry about old storage
    void resizeAndReset(int64_t newSize) {
        if (newSize <= static_cast<int64_t>(capacity())) {
            resize(newSize);
            return;
        }
        std::vector<T> other;
        other.reserve(newSize * 2);
        other.resize(newSize);
        static_cast<std::vector<T>&>(*this) = std::move(other);
    }

    /// TODO: implement push_back with parallel reallocation
    /// TODO: implement emplace_back with parallel reallocation

   private:
    using std::vector<T>::resize;
};

constexpr int maxThreads = 256;

struct Mesh::WorkspaceStencilLmat2Optimized {
    SimpleArrayBuffer<std::array<int, 3>> rowPositionsUnsorted;
    SimpleArrayBuffer<std::array<int, 3>> rowPositionsSorted;
    SimpleArrayBuffer<RowBlock<12 * 12>> globalRowBlocksMerged;
    std::array<SimpleArrayBuffer<RowBlock<12>>, maxThreads> rowBlocksLocals;
};

Mesh::Mesh() : workspacePtr(new WorkspaceStencilLmat2Optimized) {};
Mesh::~Mesh() {};

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
    std::terminate();   // avoid calling this function
    (void) mat;
    (void) domain;

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

static std::vector<Triplet> multyPhaseMerge(std::vector<std::vector<Triplet>>& local_vectors) {
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

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
static void blockToRowBlocks(int i_cell, int j_cell, int k_cell, const Block_t& block, [[maybe_unused]] int Nx, int Ny,
                             int Nz, double tolerance, std::vector<RowBlock<12>>& rowBlocks) {
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

                            if (std::abs(val) > tolerance) {
                                const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                                     k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                                if (!isAddedInThisRow) [[unlikely]] {
                                    rowBlocks.emplace_back(row);
                                    isAddedInThisRow = true;
                                }
                                rowBlocks.back().push_back_value(col, val);
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename RowIdx, typename ColIdx, int DIR, typename Block_t>
static void blockToTriplets(int i_cell, int j_cell, int k_cell, const Block_t& block, std::vector<Triplet>& trips,
                            [[maybe_unused]] int Nx, int Ny, int Nz, double tolerance) {
    auto vind = [&](int i, int j, int k, int d) { return d + 3 * (i * Ny * Nz + j * Nz + k); };
    for (int x1 = 0; x1 < RowIdx::size_x; ++x1)
        for (int y1 = 0; y1 < RowIdx::size_y; ++y1)
            for (int z1 = 0; z1 < RowIdx::size_z; ++z1) {
                const int row = vind(i_cell + x1 + RowIdx::offset_x, j_cell + y1 + RowIdx::offset_y,
                                     k_cell + z1 + RowIdx::offset_z, RowIdx::dir);

                for (int x2 = 0; x2 < ColIdx::size_x; ++x2)
                    for (int y2 = 0; y2 < ColIdx::size_y; ++y2)
                        for (int z2 = 0; z2 < ColIdx::size_z; ++z2) {
                            const double val = block(RowIdx::calculate(x1, y1, z1), ColIdx::calculate(x2, y2, z2), DIR);

                            if (std::abs(val) > tolerance) {
                                const int col = vind(i_cell + x2 + ColIdx::offset_x, j_cell + y2 + ColIdx::offset_y,
                                                     k_cell + z2 + ColIdx::offset_z, ColIdx::dir);
                                trips.emplace_back(Triplet{row, col, val});
                            }
                        }
            }
}

void Mesh::stencil_Lmat2(Operator& mat, const Domain& domain,
                         std::unique_ptr<WorkspaceStencilLmat2Optimized>& workspacePtr) const {
    RECORD_TIMER;
    static const int checkPeriodicity = envOptions::validationPeriodicity();
    static int counter = 0;
    const bool doCheck = counter % checkPeriodicity == 0;
    counter += 1;
    Operator ref;
    if (doCheck) {
        ref = mat;
        stencil_Lmat2_Reference(ref, domain);
    }
    stencil_Lmat2_Optimized(mat, domain, workspacePtr);

    if (doCheck) {
        timer::commonTimer timer("Check optimized algorithm");
        const bool isSameShape = checkMatrixPortraitCoincidence(mat, ref);
        if (!isSameShape) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line()
                      << "Critical error: optimized and reference algorithms provided matrices with portraits "
                      << std::endl;
        }

        const Operator diff = mat - ref;
        const double normalizedErr = diff.norm() / ref.norm();

        if (normalizedErr >= 1e-16) {
            const std::source_location location = std::source_location::current();
            std::cerr << location.file_name() << ":" << location.line()
                      << " Error between optimized and reference algorithms is too large: normalized error = "
                      << normalizedErr << " >= " << 1e-16 << std::endl;
        }
    }
}

void Mesh::stencil_Lmat2_Optimized(Operator& mat, const Domain& domain,
                                   std::unique_ptr<WorkspaceStencilLmat2Optimized>& workspacePtr) const {
    RECORD_TIMER;

    WorkspaceStencilLmat2Optimized& workspace = *workspacePtr;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    const auto size = domain.size();
    const int max_i = size.x() - 1;
    const int max_j = size.y() - 1;
    const int max_k = size.z() - 1;

    const int rows = mat.rows();
    const int nthr = std::min(maxThreads, omp_get_max_threads());

    std::array<SimpleArrayBuffer<RowBlock<12>>, maxThreads>& rowBlocksLocals = workspace.rowBlocksLocals;

    timer::commonTimer timerUnpackingNew("unpacking new");
#pragma omp parallel num_threads(nthr)
    {
        SimpleArrayBuffer<RowBlock<12>>& localRowBlocks = rowBlocksLocals[omp_get_thread_num()];
        localRowBlocks.resizeAndReset(0);
        localRowBlocks.reserve(1024 * 1024 * 4);
        timer::commonTimer timerOmp("OMP loop");
#pragma omp for collapse(3) schedule(dynamic, 8)
        for (int i = BORDER; i < max_i; ++i) {
            for (int j = BORDER; j < max_j; ++j) {
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = LmatX2[sind(i, j, k)];

                    // X component
                    blockToRowBlocks<XIndexer, XIndexer, 0>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);
                    blockToRowBlocks<XIndexer, YIndexer, 1>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);
                    blockToRowBlocks<XIndexer, ZIndexer, 2>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);

                    // Y component
                    blockToRowBlocks<YIndexer, XIndexer, 3>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);
                    blockToRowBlocks<YIndexer, YIndexer, 4>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);
                    blockToRowBlocks<YIndexer, ZIndexer, 5>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);

                    // Z component
                    blockToRowBlocks<ZIndexer, XIndexer, 6>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);
                    blockToRowBlocks<ZIndexer, YIndexer, 7>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);
                    blockToRowBlocks<ZIndexer, ZIndexer, 8>(i, j, k, block, xSize, ySize, zSize, TOL, localRowBlocks);
                }
            }
        }
    }
    timerUnpackingNew.finish();

    timer::commonTimer timerPseudoSort("pseudo sort");

    // shall be int64_t

    std::vector<int> offsets(nthr);
    offsets[0] = 0;
    int64_t unmergedRowBlocks = 0;
    for (int i = 0; i < nthr; ++i) {
        unmergedRowBlocks += std::ssize(rowBlocksLocals[i]);
    }
    for (int i = 1; i < nthr; ++i) {
        offsets[i] = std::ssize(rowBlocksLocals[i - 1]) + offsets[i - 1];
    }

    timer::commonTimer timerInitForSort("init vectors of ints");
    // Indexes of each arrays is: row, orig thread owner and orig. number in its buffer
    SimpleArrayBuffer<std::array<int, 3>>& rowPositionsUnsorted = workspace.rowPositionsUnsorted;
    SimpleArrayBuffer<std::array<int, 3>>& rowPositionsSorted = workspace.rowPositionsSorted;
    rowPositionsUnsorted.resizeAndReset(unmergedRowBlocks);
    rowPositionsSorted.resizeAndReset(unmergedRowBlocks);
    timerInitForSort.finish();

    timer::commonTimer timerSetUnsortedRowPos("set row positions unsorted");

#pragma omp parallel for num_threads(nthr)
    for (int i = 0; i < nthr; ++i) {
        timer::commonTimer timerOMP("OMP section");

        const SimpleArrayBuffer<RowBlock<12>>& localRowBlocks = rowBlocksLocals[i];

        for (int j = 0; j < std::ssize(localRowBlocks); ++j) {
            rowPositionsUnsorted[offsets[i] + j] = {localRowBlocks[j].row, i, j};
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
        const SimpleArrayBuffer<RowBlock<12>>& localRowBlocks = rowBlocksLocals[i];
        for (int j = 0; j < std::ssize(localRowBlocks); ++j) {
            std::atomic_ref<uint8_t> toUpd(nonZeroBlocks[localRowBlocks[j].row]);
            toUpd += 1;
        }
    }

    timerComputingNnzBlocks.finish();

    timer::commonTimer timerMulti("set blocks bound and zero aux array", sizeof(int) * rows * 2,
                                  timer::MeasureUnit::byte);

    std::vector<int> nonZeroBlocksOuter(rows + 1);
    nonZeroBlocksOuter[0] = 0;
    int nonZeroBlocksCount = 0;
    /// NOTE: could be parallelized
    for (int i = 0; i < rows; ++i) {
        nonZeroBlocksCount += nonZeroBlocks[i];
        nonZeroBlocks[i] = 0;
        nonZeroBlocksOuter[i + 1] = nonZeroBlocksCount;
    }

    timerMulti.finish();

    timer::commonTimer timerSetRowPosSorted(
        "set row positions sorted", sizeof(rowPositionsSorted[0]) * unmergedRowBlocks, timer::MeasureUnit::byte);
#pragma omp parallel for
    for (int i = 0; i < unmergedRowBlocks; ++i) {
        const int row = rowPositionsUnsorted[i][0];
        std::atomic_ref<uint8_t> blockIxRef(nonZeroBlocks[row]);
        const int blockIx = (blockIxRef++) + nonZeroBlocksOuter[row];
        rowPositionsSorted[blockIx] = rowPositionsUnsorted[i];
    }

    timerSetRowPosSorted.finish();

    timer::commonTimer timerSetBlocksBounds(
        "set block bounds in glob. array", sizeof(rowPositionsSorted[0]) * unmergedRowBlocks, timer::MeasureUnit::byte);

    timer::commonTimer timerAlloc("create vector");

    std::vector<std::vector<int>> blocksStartsUnmerged;
    blocksStartsUnmerged.resize(nthr);

    std::vector<int> blocksStarts;
    blocksStarts.reserve(unmergedRowBlocks);
    blocksStarts.push_back(0);
    timerAlloc.finish();

    /* Single thread reference code for the OMP sections below:
    for (int i = 0, j = 0; i != unmergedRowBlocks; i = j) {
        while (j < unmergedRowBlocks && rowPositionsSorted[i][0] == rowPositionsSorted[j][0]) {
            j += 1;
        }
        if (j != i) {
            blocksStarts.push_back(j);
        }
    }
    */

#pragma omp parallel num_threads(nthr)
    {
        timer::commonTimer ompTimer("OMP section");

        const int threads = omp_get_num_threads();
        const int tid = omp_get_thread_num();

        std::vector<int>& blocksStartsLocal = blocksStartsUnmerged[tid];
        blocksStartsLocal.reserve(unmergedRowBlocks / threads);
        if (tid == 0)
            blocksStartsLocal.push_back(0);

        int64_t start = unmergedRowBlocks * tid / threads;
        int64_t end = unmergedRowBlocks * (tid + 1) / threads;
        while (start != 0 && start < unmergedRowBlocks &&
               rowPositionsSorted[start - 1][0] == rowPositionsSorted[start][0]) {
            start += 1;
        }
        while (end != unmergedRowBlocks && end < unmergedRowBlocks &&
               rowPositionsSorted[end - 1][0] == rowPositionsSorted[end][0]) {
            end += 1;
        }

        for (int i = start, j = start; i != end; i = j) {
            while (j < end && rowPositionsSorted[i][0] == rowPositionsSorted[j][0]) {
                j += 1;
            }
            if (j != i) {
                blocksStartsLocal.push_back(j);
            }
        }
    }
    int blocksCount = 0;
    for (int i = 0; i < nthr; ++i) {
        blocksCount += std::ssize(blocksStartsUnmerged[i]);
    }
    blocksStarts.resize(blocksCount);

#pragma omp parallel num_threads(nthr)
    {
        const int tid = omp_get_thread_num();
        std::vector<int>& blocksStartsLocal = blocksStartsUnmerged[tid];
        int offset = 0;
        for (int i = 0; i < tid; ++i) {
            offset += std::ssize(blocksStartsUnmerged[i]);
        }
        std::copy(blocksStartsLocal.begin(), blocksStartsLocal.end(), blocksStarts.begin() + offset);
    }

    timerSetBlocksBounds.finish();

    timerPseudoSort.finish();

    timer::commonTimer timerMerge("merge", -1, timer::MeasureUnit::byte);

    SimpleArrayBuffer<RowBlock<12 * 12>>& globalRowBlocksMerged = workspace.globalRowBlocksMerged;
    globalRowBlocksMerged.resizeAndReset(blocksStarts.size() - 1);

    int totalNnz = 0;
    int64_t mergedNnz = 0;
#pragma omp parallel reduction(+ : totalNnz, mergedNnz)
    {
        timer::commonTimer timerOmp("OMP section");
        constexpr int maxMergedBlocks = 12 * 12;
        std::array<RowBlock<12>, maxMergedBlocks> tmpStorage;
#pragma omp for
        for (int i = 0; i < std::ssize(blocksStarts) - 1; ++i) {
            const int start = blocksStarts[i];
            const int end = blocksStarts[i + 1];
            assert(end - start <= maxMergedBlocks);
            for (int j = start; j < end; ++j) {
                const auto [row, thread, index] = rowPositionsSorted[j];
                tmpStorage[j - start] = rowBlocksLocals[thread][index];
                mergedNnz += tmpStorage[j - start].nnz;
            }
            globalRowBlocksMerged[i].mergeFromOthers(end - start, &tmpStorage[0]);
            totalNnz += globalRowBlocksMerged[i].nnz;
        }
        timerOmp.finish();
    }
    timerMerge.flat.m = (sizeof(int) + sizeof(double)) * mergedNnz;
    timerMerge.finish();

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

    const int64_t sizeOuterBytes = static_cast<int64_t>(sizeof(int)) * std::ssize(globalRowBlocksMerged);
    const int64_t sizeRestMatrix = static_cast<int64_t>(sizeof(int) + sizeof(double)) * totalNnz;
    timer::commonTimer timerFilling("filling");

    int* outer = mat.outerIndexPtr();
    int* ind = mat.innerIndexPtr();
    double* values = mat.valuePtr();

    timer::commonTimer timerFillingOuter("filling outer", sizeOuterBytes, timer::MeasureUnit::byte);
#pragma omp parallel for
    for (int i = 0; i < rows + 1; ++i) {
        outer[i] = 0;
    }

#pragma omp parallel for
    for (int i = 0; i < std::ssize(globalRowBlocksMerged); ++i) {
        outer[globalRowBlocksMerged[i].row + 1] = globalRowBlocksMerged[i].nnz;
    }

    /// NOTE: could be parallelized
    for (int i = 1; i < rows + 1; ++i) {
        outer[i] = outer[i - 1] + outer[i];
    }
    timerFillingOuter.finish();

    timer::commonTimer timerFillIndsAndVals("fill rest data", sizeRestMatrix, timer::MeasureUnit::byte);
#pragma omp parallel for
    for (int i = 0; i < std::ssize(globalRowBlocksMerged); ++i) {
        RowBlock<12 * 12>& rowBlock = globalRowBlocksMerged[i];
        const int row = rowBlock.row;
        const int start = outer[row];
        const int size = outer[row + 1] - start;

        for (int j = 0; j < size; ++j) {
            values[start + j] = rowBlock.values[j];
        }
        for (int j = 0; j < size; ++j) {
            ind[start + j] = rowBlock.columns[j];
        }
    }
}

void Mesh::stencil_Lmat2_Reference(Operator& mat, const Domain& domain) const {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    static int check_count = 0;
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

    // Локальные векторы для каждого потока с предварительным резервированием памяти.
    std::vector<std::vector<Triplet>> local_vectors(num_threads);

    timer::commonTimer timer1("parallel sections");

#pragma omp parallel num_threads(num_threads)
    {
        const int tid = omp_get_thread_num();
        // First touch - выделяем память в том же потоке, который будет её использовать
        local_vectors[tid].reserve(estimated_per_thread);
    }
#pragma omp parallel num_threads(num_threads)
    {
        timer::flatTimer timerLocal("main OMP section");

        const int tid = omp_get_thread_num();
        std::vector<Triplet>& localVector = local_vectors[tid];

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX2.non_zeros[sind(i, j, k)])
                        continue;
                    const Block& block = LmatX2[sind(i, j, k)];

                    // X component
                    blockToTriplets<XIndexer, XIndexer, 0>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<XIndexer, YIndexer, 1>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<XIndexer, ZIndexer, 2>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);

                    // Y component
                    blockToTriplets<YIndexer, XIndexer, 3>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<YIndexer, YIndexer, 4>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<YIndexer, ZIndexer, 5>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);

                    // Z component
                    blockToTriplets<ZIndexer, XIndexer, 6>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<ZIndexer, YIndexer, 7>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<ZIndexer, ZIndexer, 8>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                }

        // Сортировка локального вектора по (row, col)
        timer::flatTimer timerDuplications("remove duplications", localVector.size());
        timer::flatTimer timerSort("sort local triplet vector", localVector.size());
        std::sort(localVector.begin(), localVector.end(), compareTriplets);
        timerSort.finish();

        // Устранение дубликатов: проход по отсортированному вектору, складываем значения для одинаковых ключей
        if (!localVector.empty()) {
            size_t index = 0;
            for (size_t j = 1; j < localVector.size(); ++j) {
                if (localVector[index].row() == localVector[j].row() &&
                    localVector[index].col() == localVector[j].col()) {
                    localVector[index].value() += localVector[j].value();
                } else {
                    ++index;
                    localVector[index] = localVector[j];
                }
            }
            localVector.resize(index + 1);
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

    timer3.finish();

    optimizedSetFromSortedTriplets(mat, test);

    LOG_STEP("Matrix L (block) was created." << " trips size: " << test.size() << std::endl);
}

// TODO implement parallelMultiwayMergeSort and merge to convert block to csr
// (it is dublicated)
void Mesh::stencil_Lmat2_NGP(Operator& mat, const Domain& domain) {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;

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

    // Локальные векторы для каждого потока с предварительным резервированием памяти.
    std::vector<std::vector<Triplet>> local_vectors(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        local_vectors[t].reserve(estimated_per_thread);
    }

#pragma omp parallel num_threads(num_threads)
    {
        const int tid = omp_get_thread_num();

        std::vector<Triplet>& localVector = local_vectors[tid];

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!LmatX_NGP.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = LmatX_NGP[sind(i, j, k)];

                    // X component
                    blockToTriplets<XIndexerNGP, XIndexerNGP, 0>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<XIndexerNGP, YIndexerNGP, 1>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<XIndexerNGP, ZIndexerNGP, 2>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);

                    // Y component
                    blockToTriplets<YIndexerNGP, XIndexerNGP, 3>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<YIndexerNGP, YIndexerNGP, 4>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<YIndexerNGP, ZIndexerNGP, 5>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);

                    // Z component
                    blockToTriplets<ZIndexerNGP, XIndexerNGP, 6>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<ZIndexerNGP, YIndexerNGP, 7>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<ZIndexerNGP, ZIndexerNGP, 8>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                }

        // Сортировка локального вектора по (row, col)
        std::sort(localVector.begin(), localVector.end(), compareTriplets);

        // Устранение дубликатов: проход по отсортированному вектору, складываем значения для одинаковых ключей
        if (!localVector.empty()) {
            size_t index = 0;
            for (size_t j = 1; j < localVector.size(); ++j) {
                if (localVector[index].row() == localVector[j].row() &&
                    localVector[index].col() == localVector[j].col()) {
                    localVector[index].value() += localVector[j].value();
                } else {
                    ++index;
                    localVector[index] = localVector[j];
                }
            }
            localVector.resize(index + 1);
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

    // В non_empty[0] теперь находится глобальный вектор, уже отсортированный и с устранёнными дубликатами.
    const std::vector<Triplet>& trips = non_empty.empty() ? std::vector<Triplet>() : non_empty.front();

    optimizedSetFromSortedTriplets(mat, trips);
    LOG_STEP("Matrix L (block) was created." << " trips size: " << trips.size() << std::endl);
}

void Mesh::stencil_curlB(Operator& mat, const Domain& domain, BoundaryConditionHandler& bc_handler) {
    RECORD_TIMER;

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
    RECORD_TIMER;

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
    RECORD_TIMER;

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
    RECORD_TIMER;

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

template <typename IndexerX, typename IndexerY, typename IndexerZ, typename MatrixType>
void Mesh::convert_block_to_crs_format(MatrixType bmatrix, Operator& mat, const Domain& domain) {
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;

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

    // Локальные векторы для каждого потока с предварительным резервированием памяти.
    std::vector<std::vector<Triplet>> local_vectors(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        local_vectors[t].reserve(estimated_per_thread);
    }
#pragma omp parallel num_threads(num_threads)
    {
        const int tid = omp_get_thread_num();

        std::vector<Triplet>& localVector = local_vectors[tid];

#pragma omp for collapse(3) schedule(dynamic, 8) nowait
        for (int i = BORDER; i < max_i; ++i)
            for (int j = BORDER; j < max_j; ++j)
                for (int k = BORDER; k < max_k; ++k) {
                    if (!bmatrix.non_zeros[sind(i, j, k)])
                        continue;
                    const auto& block = bmatrix[sind(i, j, k)];

                    // X component
                    blockToTriplets<IndexerX, IndexerX, 0>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<IndexerX, IndexerY, 1>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<IndexerX, IndexerZ, 2>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);

                    // Y component
                    blockToTriplets<IndexerY, IndexerX, 3>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<IndexerY, IndexerY, 4>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<IndexerY, IndexerZ, 5>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);

                    // Z component
                    blockToTriplets<IndexerZ, IndexerX, 6>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<IndexerZ, IndexerY, 7>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                    blockToTriplets<IndexerZ, IndexerZ, 8>(i, j, k, block, localVector, xSize, ySize, zSize, TOL);
                }

        // Сортировка локального вектора по (row, col)
        std::sort(localVector.begin(), localVector.end(), compareTriplets);

        // Устранение дубликатов: проход по отсортированному вектору, складываем значения для одинаковых ключей
        if (!localVector.empty()) {
            size_t index = 0;
            for (size_t j = 1; j < localVector.size(); ++j) {
                if (localVector[index].row() == localVector[j].row() &&
                    localVector[index].col() == localVector[j].col()) {
                    localVector[index].value() += localVector[j].value();
                } else {
                    ++index;
                    localVector[index] = localVector[j];
                }
            }
            localVector.resize(index + 1);
        }
    }

    // Собираем только непустые локальные векторы для дальнейшего слияния
    std::vector<std::vector<Triplet>> non_empty;
    for (auto& v : local_vectors) {
        if (!v.empty()) {
            non_empty.push_back(std::move(v));
        }
    }

    // Многофазное слияние: объединяем пары отсортированных векторов параллельно с предварительным резервированием
    // памяти
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

    // В non_empty[0] теперь находится глобальный вектор, уже отсортированный и с устранёнными дубликатами.
    const std::vector<Triplet>& trips = non_empty.empty() ? std::vector<Triplet>() : non_empty.front();

    // mat.setFromTriplets(trips.begin(), trips.end());
    optimizedSetFromSortedTriplets(mat, trips);
}

template void Mesh::convert_block_to_crs_format<XIndexerNGP, YIndexerNGP, ZIndexerNGP>(BlockMatrixNGP bmatrix,
                                                                                       Operator& mat,
                                                                                       const Domain& domain);

template void Mesh::convert_block_to_crs_format<XIndexer, YIndexer, ZIndexer>(BlockMatrix bmatrix, Operator& mat,
                                                                              const Domain& domain);
