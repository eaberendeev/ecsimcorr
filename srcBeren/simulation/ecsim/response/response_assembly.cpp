#include <atomic>

#include "log_macros.h"
#include "memory.h"
#include "pmms.hpp"
#include "response.h"
#include "timer.h"

namespace ecsim_response {
constexpr int maxThreads = 256;
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

/* reference single-thread code for the whole function:
for (int i = 0; i < rows; ++i) {
    nonZeroBlocksCount += nonZeroBlocks[i];
    nonZeroBlocks[i] = 0;
    nonZeroBlocksOuter[i + 1] = nonZeroBlocksCount;
}
*/
static void setBlockBoundsAndZeroAuxArray(int rows, SmartPtr<uint8_t>& nonZeroBlocks,
                                          SmartPtr<int>& nonZeroBlocksOuter) {
    RECORD_TIMER_PARAMS(sizeof(int) * rows * 2, timer::MeasureUnit::byte);

    nonZeroBlocksOuter[0] = 0;
    const int nthr = omp_get_max_threads();
    std::vector<int> partialSums(nthr);
#pragma omp parallel num_threads(nthr)
    {
        checkNumThreads(nthr);

        const int tid = omp_get_thread_num();
        const auto [start, end] = arrayDivision(rows, tid, nthr);

        int partialSum = 0;
        for (int i = start; i < end; ++i) {
            partialSum += nonZeroBlocks[i];
            nonZeroBlocks[i] = 0;
            nonZeroBlocksOuter[i + 1] = partialSum;
        }
        partialSums[tid] = partialSum;
    }

    for (int i = 1; i < nthr; ++i) {
        partialSums[i] += partialSums[i - 1];
    }

#pragma omp parallel num_threads(nthr)
    {
        checkNumThreads(nthr);
        const int tid = omp_get_thread_num();
        if (tid != 0) {
            const auto [start, end] = arrayDivision(rows, tid, nthr);
            const int partialSum = partialSums[tid - 1];

            for (int i = start; i < end; ++i) {
                nonZeroBlocksOuter[i + 1] += partialSum;
            }
        }
    }
}

/* Reference single-thread code:
for (int i = 0; i < rows; ++i) {
    if (nonZeroBlocksOuter[i + 1] != nonZeroBlocksOuter[i]) {
        blocksStartsTest2.push_back(nonZeroBlocksOuter[i + 1]);
    }
}
*/
static void setBlocksBounds(int rows, SmartPtr<int>& nonZeroBlocksOuter, std::vector<int>& blocksBounds) {
    RECORD_TIMER_PARAMS(sizeof(int) * (rows + 1), timer::MeasureUnit::byte);
    assert(nonZeroBlocksOuter.size == rows + 1);

    const int nthr = omp_get_max_threads();

    std::vector<std::vector<int>> blockStartsDistr(nthr);

    int blockCount = 0;

#pragma omp parallel num_threads(nthr) reduction(+ : blockCount)
    {
        checkNumThreads(nthr);

        const int64_t tid = omp_get_thread_num();
        const auto [start, end] = arrayDivision(rows, tid, nthr);

        std::vector<int>& toUpd = blockStartsDistr[tid];
        toUpd.reserve(end - start);

        if (tid == 0) {
            toUpd.push_back(0);
        }

        for (int j = start; j < end; ++j) {
            if (nonZeroBlocksOuter[j + 1] != nonZeroBlocksOuter[j]) {
                toUpd.push_back(nonZeroBlocksOuter[j + 1]);
            }
        }

        blockCount = std::ssize(toUpd);
    }

    timer::commonTimer timerBlockStartsResize("blocksBounds.resize(blockCount)");
    blocksBounds.resize(blockCount);
    timerBlockStartsResize.finish();

#pragma omp parallel num_threads(nthr)
    {
        checkNumThreads(nthr);
        const int tid = omp_get_thread_num();
        int offset = 0;
        for (int j = 0; j < tid; ++j) {
            offset += std::ssize(blockStartsDistr[j]);
        }

        std::copy(blockStartsDistr[tid].begin(), blockStartsDistr[tid].end(), blocksBounds.begin() + offset);
    }
}

/* Reference single-thread code:
    for (int i = 0; i < rows + 1; ++i) {
        outer[i] = 0;
    }
    for (int i = 0; i < std::ssize(globalRowBlocksMerged); ++i) {
        outer[globalRowBlocksMerged[i].row + 1] = globalRowBlocksMerged[i].nnz;
    }
    for (int i = 0; i < rows; ++i) {
        outer[i + 1] = outer[i + 1] + outer[i];
    }
*/
static void fillOuter(const int rows, int* outer, const SimpleArrayBuffer<RowBlock<12 * 12>>& globalRowBlocksMerged) {
    const int64_t sizeOuterBytes = sizeof(int) * std::ssize(globalRowBlocksMerged);
    RECORD_TIMER_PARAMS(sizeOuterBytes, timer::MeasureUnit::byte);

#pragma omp parallel for
    for (int i = 0; i < rows + 1; ++i) {
        outer[i] = 0;
    }

#pragma omp parallel for
    for (int i = 0; i < std::ssize(globalRowBlocksMerged); ++i) {
        outer[globalRowBlocksMerged[i].row + 1] = globalRowBlocksMerged[i].nnz;
    }

    const int nthr = omp_get_max_threads();
    if (rows < nthr * 2) {
        // reference single-thread code for the parallel code below:
        for (int i = 0; i < rows; ++i) {
            outer[i + 1] = outer[i + 1] + outer[i];
        }
        return;
    }

#pragma omp parallel num_threads(nthr)
    {
        checkNumThreads(nthr);
        const int tid = omp_get_thread_num();
        const auto [start, end] = arrayDivision(rows + 1, tid, nthr);
        for (int i = start; i < end - 1; ++i) {
            outer[i + 1] = outer[i + 1] + outer[i];
        }
    }

    std::vector<int> partialSums(nthr);
    for (int tid = 0; tid < nthr; ++tid) {
        const auto [start, end] = arrayDivision(rows + 1, tid, nthr);
        (void) start;
        partialSums[tid] = outer[end - 1];
    }

    for (int i = 1; i < nthr; ++i) {
        partialSums[i] += partialSums[i - 1];
    }

#pragma omp parallel num_threads(nthr)
    {
        checkNumThreads(nthr);
        const int tid = omp_get_thread_num();
        if (tid != 0) {
            const auto [start, end] = arrayDivision(rows + 1, tid, nthr);
            const int partialSum = partialSums[tid - 1];
            for (int i = start; i < end; ++i) {
                outer[i] += partialSum;
            }
        }
    }
}

void assemble_matrix(Operator& mat, const RowBlocks& row_blocks, Workspace& workspace) {
    RECORD_TIMER;
    const RowBlocks& rowBlocksLocals = row_blocks;
    assert(std::ssize(rowBlocksLocals) == omp_get_max_threads());

    const int rows = mat.rows();
    const int nthr = omp_get_max_threads();

    timer::commonTimer timerPseudoSort("pseudo sort");

    SmartPtr<int> offsets(nthr);
    offsets[0] = 0;
    int64_t unmergedRowBlocks = 0;
    for (int i = 0; i < nthr; ++i) {
        unmergedRowBlocks += std::ssize(rowBlocksLocals[i]);
    }
    for (int i = 1; i < nthr; ++i) {
        offsets[i] = std::ssize(rowBlocksLocals[i - 1]) + offsets[i - 1];
    }

    timer::commonTimer timerInitForSort("Reallocate int workspaces");
    // Indexes of each arrays is: row, orig thread owner and orig. number in its buffer
    SimpleArrayBuffer<std::array<int, 3>>& rowPositionsUnsorted = workspace.row_positions_unsorted;
    SimpleArrayBuffer<std::array<int, 3>>& rowPositionsSorted = workspace.row_positions_sorted;
    rowPositionsUnsorted.resizeAndReset(unmergedRowBlocks);
    rowPositionsSorted.resizeAndReset(unmergedRowBlocks);
    timerInitForSort.finish();

    timer::commonTimer timerSetUnsortedRowPos(
        "set row positions unsorted", sizeof(rowPositionsUnsorted[0]) * unmergedRowBlocks, timer::MeasureUnit::byte);

#pragma omp parallel for num_threads(nthr)
    for (int i = 0; i < nthr; ++i) {
        timer::commonTimer timerOMP("OMP section");

        const SimpleArrayBuffer<RowBlock<36>>& localRowBlocks = rowBlocksLocals[i];

        for (int j = 0; j < std::ssize(localRowBlocks); ++j) {
            rowPositionsUnsorted[offsets[i] + j] = {localRowBlocks[j].row, i, j};
        }
    }

    timerSetUnsortedRowPos.finish();

    timer::commonTimer timerComputingNnzBlocks("compute nnz block count in rows", (sizeof(int) + 1) * unmergedRowBlocks,
                                               timer::MeasureUnit::byte);

    SmartPtr<uint8_t> nonZeroBlocks(rows);
#pragma omp parallel for
    for (int i = 0; i < rows; ++i) {
        nonZeroBlocks[i] = 0;
    }

#pragma omp parallel
    {
        timer::flatTimer timerOMP("OMP Section");
#pragma omp for
        for (int i = 0; i < std::ssize(rowPositionsUnsorted); ++i) {
            std::atomic_ref<uint8_t> toUpd(nonZeroBlocks[rowPositionsUnsorted[i][0]]);
            const uint8_t oldVal = toUpd.fetch_add(1);
            assert(oldVal < 255);
            (void) oldVal;
        }
    }

    timerComputingNnzBlocks.finish();

    timer::commonTimer timerMulti("set blocks bound and zero aux array", sizeof(int) * rows * 2,
                                  timer::MeasureUnit::byte);

    SmartPtr<int> nonZeroBlocksOuter(rows + 1);
    setBlockBoundsAndZeroAuxArray(rows, nonZeroBlocks, nonZeroBlocksOuter);

    timerMulti.finish();

    timer::commonTimer timerSetRowPosSorted(
        "set row positions sorted", sizeof(rowPositionsSorted[0]) * unmergedRowBlocks, timer::MeasureUnit::byte);
#pragma omp parallel
    {
        timer::flatTimer timerOMP("OMP Section");
#pragma omp for
        for (int i = 0; i < unmergedRowBlocks; ++i) {
            const int row = rowPositionsUnsorted[i][0];
            std::atomic_ref<uint8_t> blockIxRef(nonZeroBlocks[row]);
            const int blockIx = (blockIxRef++) + nonZeroBlocksOuter[row];
            rowPositionsSorted[blockIx] = rowPositionsUnsorted[i];
        }
    }

    timerSetRowPosSorted.finish();

    // blocksBounds[i] and [i+1] are bounds of -th sequential sorted row-blocks for the same row
    std::vector<int> blocksBounds;
    setBlocksBounds(rows, nonZeroBlocksOuter, blocksBounds);

    timerPseudoSort.finish();

    timer::commonTimer timerMerge("merge", -1, timer::MeasureUnit::byte);

    SimpleArrayBuffer<RowBlock<12 * 12>>& globalRowBlocksMerged = workspace.merged_rows;
    globalRowBlocksMerged.resizeAndReset(blocksBounds.size() - 1);

    int totalNnz = 0;
    int64_t mergedNnz = 0;
#pragma omp parallel reduction(+ : totalNnz, mergedNnz)
    {
        timer::commonTimer timerOmp("OMP section");
#pragma omp for
        for (int i = 0; i < std::ssize(blocksBounds) - 1; ++i) {
            const int start = blocksBounds[i];
            const int end = blocksBounds[i + 1];
            assert(start <= end);
            if (end == start) {
                continue;
            }
            RowBlock<36> tmpStorage[end - start];
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
    fillOuter(rows, outer, globalRowBlocksMerged);
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

void assemble_reference(Operator& mat, const Domain& domain, const BlockMatrix& reference_blocks) {
    const BlockMatrix& LmatX2 = reference_blocks;
    const int xSize = domain.size().x(), ySize = domain.size().y(), zSize = domain.size().z();
    RECORD_TIMER;

    constexpr double TOL = 1e-16;
    constexpr int BORDER = 1;
    static int check_count = 0;
    const Vector3I size = domain.size();
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
                    if (!LmatX2.non_zeros[domain.sind(i, j, k)])
                        continue;
                    const Block& block = LmatX2[domain.sind(i, j, k)];

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

}   // namespace ecsim_response
