
#pragma once

#include <cassert>
#include <sstream>

#include "Eigen/Sparse"
#include "timer.h"
#include "types.h"
#include "util.h"

namespace thread_partitioned_mat_impl {
/* Used to find best partitions of thread-partitioned matrix
 */
template <typename T>
inline int findBestPos(const Eigen::SparseMatrix<T, MAJOR>& A, int blockId, int numBlocks) {
    RECORD_TIMER;
    if (blockId == 0) {
        return 0;
    }
    // this check is crucial, since provides correct index of the last block, for case when last rows of matrix
    // A are empty
    if (blockId == numBlocks) {
        return A.rows();
    }

    const int nnz = A.nonZeros();
    const int bestPos = static_cast<int64_t>(nnz) * blockId / numBlocks;

    const int* outer = A.outerIndexPtr();

    for (int i = 0; i < A.rows(); ++i) {
        if (outer[i] <= bestPos && bestPos <= outer[i + 1]) {
            return i;
        }
    }
    return std::numeric_limits<int>::min();
}

inline int mergeSorted(const std::vector<std::array<int, 256>>& inArrays, std::vector<int> sizes,
                       std::array<int, 256>& res) {
    RECORD_TIMER;
    std::vector<int> offsets(sizes.size(), 0);

    int pos = 0;
    while (true) {
        int minVal = std::numeric_limits<int>::max();
        for (int i = 0; i < std::ssize(inArrays); ++i) {
            if (offsets[i] < sizes[i]) {
                minVal = std::min(minVal, inArrays[i][offsets[i]]);
            }
        }

        if (minVal == std::numeric_limits<int>::max()) {
            break;
        }

        assert(pos < 256);
        res[pos] = minVal;
        pos += 1;

        for (int i = 0; i < std::ssize(inArrays); ++i) {
            if (offsets[i] < sizes[i] && minVal == inArrays[i][offsets[i]]) {
                offsets[i] += 1;
            }
        }
    }

    return pos;
}

/* Used for greedy matrix representation
 */
inline int SetOffsets(const Eigen::SparseMatrix<double, MAJOR>& A, std::array<int, 256>& colOffsets) {
    constexpr int64_t sizeofElem = sizeof(A.innerIndexPtr()[0]);
    RECORD_TIMER_PARAMS(A.nonZeros() * sizeofElem, timer::MeasureUnit::byte);

    const int rows = A.rows();
    const int* inner = A.innerIndexPtr();
    const int* outer = A.outerIndexPtr();

    std::vector<std::array<int, 256>> localOffsets(omp_get_max_threads());
    std::vector<int> offsetsCountGlobal(omp_get_max_threads());

#pragma omp parallel
    {
        timer::commonTimer timerOmp("OMP section 1");
        int offsetsCount = 0;
        std::array<int, 256>& offsets = localOffsets[omp_get_thread_num()];

#pragma omp for schedule(dynamic, 16 * 1024)
        for (int i = 0; i < rows; ++i) {
            for (int j = outer[i]; j < outer[i + 1]; ++j) {
                const int diff = inner[j] - i;
                if (std::find(offsets.begin(), offsets.begin() + offsetsCount, diff) ==
                    offsets.begin() + offsetsCount) {
                    if (offsetsCount >= 256) {
                        throw std::runtime_error("can't handle more than 256 offsets due to byte limit");
                    }
                    offsets[offsetsCount] = diff;
                    offsetsCount += 1;
                }
            }
        }
        std::sort(offsets.begin(), offsets.begin() + offsetsCount);
#pragma omp atomic write
        offsetsCountGlobal[omp_get_thread_num()] = offsetsCount;
    }

    const int offsetsCount = thread_partitioned_mat_impl::mergeSorted(localOffsets, offsetsCountGlobal, colOffsets);
    return offsetsCount;
}

}   // namespace thread_partitioned_mat_impl

/* Utilizes numa-aware storage for sparse matrices: we suppose, what newly allocated memory when in touched goes into
 * physical memory, local to writing thread. So, we re-store original sparse matrix in such way, that in smpv each
 * thread will access only data, which is located at his num node
 */
template <typename T>
struct ThreadPartitionedSparseMatrix {
   private:
    struct SparseSubMatrix;

   public:
    template <typename other_t>
    ThreadPartitionedSparseMatrix(const Eigen::SparseMatrix<other_t, MAJOR>& A)
        : nthr(omp_get_max_threads()), nnz(A.nonZeros()), matrices(nthr) {
        constexpr int64_t sizeofElem =
            (std::max(sizeof(T), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));
        RECORD_TIMER_PARAMS(A.nonZeros() * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel num_threads(nthr)
        {
            if (omp_get_num_threads() != nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }
            matrices[omp_get_thread_num()].init(A);
        }
    }

    template <typename VectorType>
    friend void spmv(const ThreadPartitionedSparseMatrix& A, const VectorType& v, VectorType& res) {
        constexpr int64_t sizeofElem = (sizeof(T) + sizeof(A.matrices[0].innerIndexes[0]));
        RECORD_TIMER_PARAMS(A.nnz * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel num_threads(A.nthr)
        {
            if (omp_get_num_threads() != A.nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(A.nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }

            const typename ThreadPartitionedSparseMatrix<T>::SparseSubMatrix& localMat =
                A.matrices[omp_get_thread_num()];

            timer::flatTimer timerOMP("OMP section",
                                      (localMat.outerIndexes.back() - localMat.outerIndexes.front()) * sizeofElem,
                                      timer::MeasureUnit::byte);

            for (int i = localMat.rowStart; i < localMat.rowEnd; ++i) {
                double sum = 0.0;
                const int localRow = i - localMat.rowStart;
                const int dataOffset = localMat.outerIndexes[0];
                const int start = localMat.outerIndexes[localRow] - dataOffset;
                const int end = localMat.outerIndexes[localRow + 1] - dataOffset;
#pragma omp simd
                for (int j = start; j < end; ++j) {
                    sum += static_cast<double>(localMat.data[j]) * static_cast<double>(v[localMat.innerIndexes[j]]);
                }
                res[i] = static_cast<T>(sum);
            }
        }
    }

    const int nthr;
    const int nnz;   // for timings only
    std::vector<SparseSubMatrix> matrices;

   private:
    struct SparseSubMatrix {
        template <typename other_t>
        void init(const Eigen::SparseMatrix<other_t, MAJOR>& A) {
            RECORD_TIMER;
            const int numThreads = omp_get_num_threads();
            const int tid = omp_get_thread_num();

            threadOwner = tid;
            rowsGlobal = A.rows();

            const other_t* val = A.valuePtr();
            const int* inner = A.innerIndexPtr();
            const int* outer = A.outerIndexPtr();

            rowStart = thread_partitioned_mat_impl::findBestPos(A, tid, numThreads);
            rowEnd = thread_partitioned_mat_impl::findBestPos(A, tid + 1, numThreads);

            outerIndexes.resize(rowEnd - rowStart + 1);
            innerIndexes.resize(outer[rowEnd] - outer[rowStart]);
            data.resize(outer[rowEnd] - outer[rowStart]);

            for (int i = rowStart; i < rowEnd + 1; ++i) {
                outerIndexes[i - rowStart] = outer[i];
            }

            constexpr int64_t sizeofElem =
                (std::max(sizeof(T), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));
            timer::commonTimer timerCopying("copy row-block", (outer[rowEnd] - outer[rowStart]) * sizeofElem,
                                            timer::MeasureUnit::byte);
            const int dataOffset = outer[rowStart];
            for (int i = rowStart; i < rowEnd; ++i) {
                for (int j = outer[i]; j < outer[i + 1]; ++j) {
                    data[j - dataOffset] = static_cast<T>(val[j]);
                    innerIndexes[j - dataOffset] = inner[j];
                }
            }
        }

        int rowsGlobal;
        int rowStart;
        int rowEnd;
        int threadOwner;

        std::vector<int> innerIndexes;
        std::vector<int> outerIndexes;
        std::vector<T> data;
    };
};

/* The same like ThreadPartitionedSparseMatrix, but uses lower integer types for the inner indexes
 */
template <typename T>
struct GreedyThreadPartitionedSparseMatrix {
   private:
    struct GreedySparseSubMatrix;

   public:
    template <typename other_t>
    GreedyThreadPartitionedSparseMatrix(const Eigen::SparseMatrix<other_t, MAJOR>& A)
        : nthr(omp_get_max_threads()), nnz(A.nonZeros()), matrices(nthr) {
        constexpr int64_t sizeofElem =
            (std::max(sizeof(T), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));
        RECORD_TIMER_PARAMS(A.nonZeros() * sizeofElem, timer::MeasureUnit::byte);

        const int offsetsCount = thread_partitioned_mat_impl::SetOffsets(A, colOffsets);

#pragma omp parallel num_threads(nthr)
        {
            if (omp_get_num_threads() != nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }

            timer::commonTimer timerOmp("OMP section 2");
            matrices[omp_get_thread_num()].init(A, colOffsets, offsetsCount);
        }
    }

    template <typename VectorType>
    friend void spmv(const GreedyThreadPartitionedSparseMatrix& A, const VectorType& v, VectorType& res) {
        constexpr int64_t sizeofElem = (sizeof(T) + sizeof(A.matrices[0].offsetInnerIndexes[0]));
        RECORD_TIMER_PARAMS(A.nnz * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel num_threads(A.nthr)
        {
            if (omp_get_num_threads() != A.nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(A.nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }

            const typename GreedyThreadPartitionedSparseMatrix<T>::GreedySparseSubMatrix& localMat =
                A.matrices[omp_get_thread_num()];

            timer::flatTimer timerOMP("OMP section",
                                      (localMat.outerIndexes.back() - localMat.outerIndexes.front()) * sizeofElem,
                                      timer::MeasureUnit::byte);

            for (int i = localMat.rowStart; i < localMat.rowEnd; ++i) {
                double sum = 0.0;
                const int localRow = i - localMat.rowStart;
                const int dataOffset = localMat.outerIndexes[0];
                const int start = localMat.outerIndexes[localRow] - dataOffset;
                const int end = localMat.outerIndexes[localRow + 1] - dataOffset;
#pragma omp simd
                for (int j = start; j < end; ++j) {
                    const int offset = A.colOffsets[localMat.offsetInnerIndexes[j]];
                    const int col = i + offset;
                    sum += static_cast<double>(localMat.data[j]) * static_cast<double>(v[col]);
                }
                res[i] = static_cast<T>(sum);
            }
        }
    }

    std::array<int, 256> colOffsets;
    int offsetsCount;
    const int nthr;
    const int nnz;   // for timings only
    std::vector<GreedySparseSubMatrix> matrices;

   private:
    struct GreedySparseSubMatrix {
        template <typename other_t>
        void init(const Eigen::SparseMatrix<other_t, MAJOR>& A, const std::array<int, 256>& offsets, int offsetsCount) {
            RECORD_TIMER;
            const int numThreads = omp_get_num_threads();
            const int tid = omp_get_thread_num();

            threadOwner = tid;
            rowsGlobal = A.rows();

            const other_t* val = A.valuePtr();
            const int* inner = A.innerIndexPtr();
            const int* outer = A.outerIndexPtr();

            rowStart = thread_partitioned_mat_impl::findBestPos(A, tid, numThreads);
            rowEnd = thread_partitioned_mat_impl::findBestPos(A, tid + 1, numThreads);

            outerIndexes.resize(rowEnd - rowStart + 1);
            offsetInnerIndexes.resize(outer[rowEnd] - outer[rowStart]);
            data.resize(outer[rowEnd] - outer[rowStart]);

            for (int i = rowStart; i < rowEnd + 1; ++i) {
                outerIndexes[i - rowStart] = outer[i];
            }

            constexpr int64_t sizeofElem =
                (std::max(sizeof(T), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));
            timer::commonTimer timerCopying("copy row-block", (outer[rowEnd] - outer[rowStart]) * sizeofElem,
                                            timer::MeasureUnit::byte);
            const int dataOffset = outer[rowStart];
            for (int i = rowStart; i < rowEnd; ++i) {
                for (int j = outer[i]; j < outer[i + 1]; ++j) {
                    data[j - dataOffset] = static_cast<T>(val[j]);
                    const int offset = inner[j] - i;
                    const int* pos = std::find(offsets.data(), offsets.data() + offsetsCount, offset);
                    const int offsetIndex = (pos - offsets.data());
                    assert(offsetIndex < offsetsCount);
                    offsetInnerIndexes[j - dataOffset] = offsetIndex;
                }
            }
        }

        int rowsGlobal;
        int rowStart;
        int rowEnd;
        int threadOwner;

        std::vector<uint8_t> offsetInnerIndexes;
        std::vector<int> outerIndexes;
        std::vector<T> data;
    };
};

template <typename T>
struct ThreadPartitionedSparseMatrixView {
    ThreadPartitionedSparseMatrixView(const std::vector<int>& rowStartsIn, const std::vector<int>& rowEndsIn,
                                      const std::vector<std::vector<int>>& innerIndexesGlobIn,
                                      const std::vector<std::vector<int>>& outerIndexesGlobIn,
                                      const std::vector<std::vector<T>>& dataGlobIn, int nnzIn, int nthrIn)
        : rowStarts(rowStartsIn),
          rowEnds(rowEndsIn),
          innerIndexesGlob(innerIndexesGlobIn),
          outerIndexesGlob(outerIndexesGlobIn),
          dataGlob(dataGlobIn),
          nnz(nnzIn),
          nthr(nthrIn) {
    }

    template <typename VectorType>
    friend void spmv(const ThreadPartitionedSparseMatrixView& A, const VectorType& v, VectorType& res) {
        constexpr int64_t sizeofElem = (sizeof(T) + sizeof(A.innerIndexesGlob[0][0]));
        RECORD_TIMER_PARAMS(A.nnz * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel num_threads(A.nthr)
        {
            if (omp_get_num_threads() != A.nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(A.nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }

            const int tid = omp_get_thread_num();
            const std::vector<int>& innerIndexes = A.innerIndexesGlob[tid];
            const std::vector<int>& outerIndexes = A.outerIndexesGlob[tid];
            const std::vector<T>& data = A.dataGlob[tid];
            const int rowStart = A.rowStarts[tid];
            const int rowEnd = A.rowEnds[tid];

            timer::flatTimer timerOMP("OMP section", (outerIndexes.back() - outerIndexes.front()) * sizeofElem,
                                      timer::MeasureUnit::byte);

            for (int i = rowStart; i < rowEnd; ++i) {
                double sum = 0.0;
                const int localRow = i - rowStart;
                const int dataOffset = outerIndexes[0];
                const int start = outerIndexes[localRow] - dataOffset;
                const int end = outerIndexes[localRow + 1] - dataOffset;
#pragma omp simd
                for (int j = start; j < end; ++j) {
                    sum += static_cast<double>(data[j]) * static_cast<double>(v[innerIndexes[j]]);
                }
                res[i] = static_cast<T>(sum);
            }
        }
    }

    const std::vector<int>& rowStarts;
    const std::vector<int>& rowEnds;
    const std::vector<std::vector<int>>& innerIndexesGlob;
    const std::vector<std::vector<int>>& outerIndexesGlob;
    const std::vector<std::vector<T>>& dataGlob;
    int nnz;
    int nthr;
};

template <typename T1, typename T2>
struct ThreadPartitionedSparseMatrixArray {
   private:
    struct SparseSubMatrix;

   public:
    template <typename other_t>
    ThreadPartitionedSparseMatrixArray(const Eigen::SparseMatrix<other_t, MAJOR>& A) {
        constexpr int64_t sizeofElem =
            (std::max(sizeof(T1) + sizeof(T2), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));

        RECORD_TIMER_PARAMS(A.nonZeros() * sizeofElem, timer::MeasureUnit::byte);
        nnz = A.nonZeros();

        nthr = omp_get_max_threads();
        rowStarts.resize(nthr);
        rowEnds.resize(nthr);
        innerIndexesGlob.resize(nthr);
        outerIndexesGlob.resize(nthr);
        dataGlob1.resize(nthr);
        dataGlob2.resize(nthr);

#pragma omp parallel num_threads(nthr)
        {
            if (omp_get_num_threads() != nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }

            int tid = omp_get_thread_num();
            std::vector<int>& innerIndexes = innerIndexesGlob[tid];
            std::vector<int>& outerIndexes = outerIndexesGlob[tid];
            std::vector<T1>& data1 = dataGlob1[tid];
            std::vector<T2>& data2 = dataGlob2[tid];

            const other_t* val = A.valuePtr();
            const int* inner = A.innerIndexPtr();
            const int* outer = A.outerIndexPtr();

            const int rowStart = thread_partitioned_mat_impl::findBestPos(A, tid, nthr);
            const int rowEnd = thread_partitioned_mat_impl::findBestPos(A, tid + 1, nthr);
            rowStarts[tid] = rowStart;
            rowEnds[tid] = rowEnd;

            outerIndexes.resize(rowEnd - rowStart + 1);
            innerIndexes.resize(outer[rowEnd] - outer[rowStart]);
            data1.resize(outer[rowEnd] - outer[rowStart]);
            data2.resize(outer[rowEnd] - outer[rowStart]);

            for (int i = rowStart; i < rowEnd + 1; ++i) {
                outerIndexes[i - rowStart] = outer[i];
            }

            timer::commonTimer timerCopying("copy row-block", (outer[rowEnd] - outer[rowStart]) * sizeofElem,
                                            timer::MeasureUnit::byte);
            const int dataOffset = outer[rowStart];
            for (int i = rowStart; i < rowEnd; ++i) {
                for (int j = outer[i]; j < outer[i + 1]; ++j) {
                    const other_t readVal = val[j];
                    data1[j - dataOffset] = static_cast<T1>(readVal);
                    data2[j - dataOffset] = static_cast<T2>(readVal);
                    innerIndexes[j - dataOffset] = inner[j];
                }
            }
        }
    }

    template <typename T>
    ThreadPartitionedSparseMatrixView<T> get() const {
        static_assert(std::is_same_v<T, T1> || std::is_same_v<T, T2>);
        if constexpr (std::is_same_v<T, T1>) {
            return ThreadPartitionedSparseMatrixView<T>(rowStarts, rowEnds, innerIndexesGlob, outerIndexesGlob,
                                                        dataGlob1, nnz, nthr);
        } else if constexpr (std::is_same_v<T, T2>)
            return ThreadPartitionedSparseMatrixView<T>(rowStarts, rowEnds, innerIndexesGlob, outerIndexesGlob,
                                                        dataGlob2, nnz, nthr);
    }

    std::vector<int> rowStarts;
    std::vector<int> rowEnds;
    std::vector<std::vector<int>> innerIndexesGlob;
    std::vector<std::vector<int>> outerIndexesGlob;
    std::vector<std::vector<T1>> dataGlob1;
    std::vector<std::vector<T2>> dataGlob2;
    int nnz;
    int nthr;
};

template <typename T>
struct GreedyThreadPartitionedSparseMatrixView {
    GreedyThreadPartitionedSparseMatrixView(const std::vector<int>& rowStartsIn, const std::vector<int>& rowEndsIn,
                                            const std::vector<std::vector<uint8_t>>& offsetInnerIndexesGlobIn,
                                            const std::vector<std::vector<int>>& outerIndexesGlobIn,
                                            const std::vector<std::vector<T>>& dataGlobIn,
                                            const std::array<int, 256>& colOffsetsIn, int nnzIn, int nthrIn)
        : rowStarts(rowStartsIn),
          rowEnds(rowEndsIn),
          offsetInnerIndexesGlob(offsetInnerIndexesGlobIn),
          outerIndexesGlob(outerIndexesGlobIn),
          dataGlob(dataGlobIn),
          colOffsets(colOffsetsIn),
          nnz(nnzIn),
          nthr(nthrIn) {
    }

    template <typename VectorType>
    friend void spmv(const GreedyThreadPartitionedSparseMatrixView& A, const VectorType& v, VectorType& res) {
        constexpr int64_t sizeofElem = (sizeof(T) + sizeof(A.offsetInnerIndexesGlob[0][0]));
        RECORD_TIMER_PARAMS(A.nnz * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel num_threads(A.nthr)
        {
            if (omp_get_num_threads() != A.nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(A.nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }

            const int tid = omp_get_thread_num();
            const std::vector<uint8_t>& offsetInnerIndexes = A.offsetInnerIndexesGlob[tid];
            const std::vector<int>& outerIndexes = A.outerIndexesGlob[tid];
            const std::vector<T>& data = A.dataGlob[tid];
            const int rowStart = A.rowStarts[tid];
            const int rowEnd = A.rowEnds[tid];

            timer::flatTimer timerOMP("OMP section", (outerIndexes.back() - outerIndexes.front()) * sizeofElem,
                                      timer::MeasureUnit::byte);

            for (int i = rowStart; i < rowEnd; ++i) {
                double sum = 0.0;
                const int localRow = i - rowStart;
                const int dataOffset = outerIndexes[0];
                const int start = outerIndexes[localRow] - dataOffset;
                const int end = outerIndexes[localRow + 1] - dataOffset;
#pragma omp simd
                for (int j = start; j < end; ++j) {
                    const int offset = A.colOffsets[offsetInnerIndexes[j]];
                    const int col = i + offset;
                    sum += static_cast<double>(data[j]) * static_cast<double>(v[col]);
                }
                res[i] = static_cast<T>(sum);
            }
        }
    }

    const std::vector<int>& rowStarts;
    const std::vector<int>& rowEnds;
    const std::vector<std::vector<uint8_t>>& offsetInnerIndexesGlob;
    const std::vector<std::vector<int>>& outerIndexesGlob;
    const std::vector<std::vector<T>>& dataGlob;
    const std::array<int, 256>& colOffsets;
    int nnz;
    int nthr;
};

template <typename T1, typename T2>
struct GreedyThreadPartitionedSparseMatrixArray {
   public:
    template <typename other_t>
    GreedyThreadPartitionedSparseMatrixArray(const Eigen::SparseMatrix<other_t, MAJOR>& A) {
        constexpr int64_t sizeofElem =
            (std::max(sizeof(T1) + sizeof(T2), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));
        RECORD_TIMER_PARAMS(A.nonZeros() * sizeofElem, timer::MeasureUnit::byte);
        nnz = A.nonZeros();

        const int offsetsCount = thread_partitioned_mat_impl::SetOffsets(A, colOffsets);
        nthr = omp_get_max_threads();
        rowStarts.resize(nthr);
        rowEnds.resize(nthr);
        offsetInnerIndexesGlob.resize(nthr);
        outerIndexesGlob.resize(nthr);
        dataGlob1.resize(nthr);
        dataGlob2.resize(nthr);

#pragma omp parallel num_threads(nthr)
        {
            if (omp_get_num_threads() != nthr) {
                throw std::runtime_error("OMP created less, rather requested: " + std::to_string(nthr) + " + " +
                                         std::to_string(omp_get_num_threads()));
            }

            int tid = omp_get_thread_num();
            std::vector<uint8_t>& offsetInnerIndexes = offsetInnerIndexesGlob[tid];
            std::vector<int>& outerIndexes = outerIndexesGlob[tid];
            std::vector<T1>& data1 = dataGlob1[tid];
            std::vector<T2>& data2 = dataGlob2[tid];

            const other_t* val = A.valuePtr();
            const int* inner = A.innerIndexPtr();
            const int* outer = A.outerIndexPtr();

            const int rowStart = thread_partitioned_mat_impl::findBestPos(A, tid, nthr);
            const int rowEnd = thread_partitioned_mat_impl::findBestPos(A, tid + 1, nthr);
            rowStarts[tid] = rowStart;
            rowEnds[tid] = rowEnd;

            outerIndexes.resize(rowEnd - rowStart + 1);
            offsetInnerIndexes.resize(outer[rowEnd] - outer[rowStart]);
            data1.resize(outer[rowEnd] - outer[rowStart]);
            data2.resize(outer[rowEnd] - outer[rowStart]);

            for (int i = rowStart; i < rowEnd + 1; ++i) {
                outerIndexes[i - rowStart] = outer[i];
            }

            timer::commonTimer timerCopying("copy row-block", (outer[rowEnd] - outer[rowStart]) * sizeofElem,
                                            timer::MeasureUnit::byte);
            const int dataOffset = outer[rowStart];
            for (int i = rowStart; i < rowEnd; ++i) {
                for (int j = outer[i]; j < outer[i + 1]; ++j) {
                    const other_t readVal = val[j];
                    data1[j - dataOffset] = static_cast<T1>(readVal);
                    data2[j - dataOffset] = static_cast<T2>(readVal);

                    const int offset = inner[j] - i;
                    const int* pos = std::find(colOffsets.data(), colOffsets.data() + offsetsCount, offset);
                    const int offsetIndex = (pos - colOffsets.data());
                    assert(offsetIndex < offsetsCount);
                    offsetInnerIndexes[j - dataOffset] = offsetIndex;
                }
            }
        }
    }

    template <typename T>
    GreedyThreadPartitionedSparseMatrixView<T> get() const {
        static_assert(std::is_same_v<T, T1> || std::is_same_v<T, T2>);
        if constexpr (std::is_same_v<T, T1>) {
            return GreedyThreadPartitionedSparseMatrixView<T>(rowStarts, rowEnds, offsetInnerIndexesGlob,
                                                              outerIndexesGlob, dataGlob1, colOffsets, nnz, nthr);
        } else if constexpr (std::is_same_v<T, T2>)
            return GreedyThreadPartitionedSparseMatrixView<T>(rowStarts, rowEnds, offsetInnerIndexesGlob,
                                                              outerIndexesGlob, dataGlob2, colOffsets, nnz, nthr);
    }

    std::array<int, 256> colOffsets;
    std::vector<int> rowStarts;
    std::vector<int> rowEnds;
    std::vector<std::vector<uint8_t>> offsetInnerIndexesGlob;
    std::vector<std::vector<int>> outerIndexesGlob;
    std::vector<std::vector<T1>> dataGlob1;
    std::vector<std::vector<T2>> dataGlob2;
    int nnz;
    int nthr;
};
