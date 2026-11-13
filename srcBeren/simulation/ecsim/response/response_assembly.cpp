#include <atomic>

#include "response.h"
#include "runtime/profiler.h"
#include "utils/memory.h"
#include "utils/util.h"

namespace ecsim_response {
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

void assemble_matrix(Operator& mat, const RowBlocks& rowBlocksLocals, Workspace& workspace) {
    RECORD_TIMER;
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

}   // namespace ecsim_response
