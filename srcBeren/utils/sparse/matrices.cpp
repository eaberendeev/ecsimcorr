#include <array>
#include <vector>

#include "containers.h"
#include "sparse.h"

template <typename T>
template <typename other_t>
ThreadPartitionedSparseMatrix<T>::ThreadPartitionedSparseMatrix(const Eigen::SparseMatrix<other_t, MAJOR>& A)
    : nthr(omp_get_max_threads()), nnz(A.nonZeros()), matrices(nthr) {
    constexpr int64_t sizeofElem =
        (std::max(sizeof(T), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));
    RECORD_TIMER_PARAMS(A.nonZeros() * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel
    {
        matrices[omp_get_thread_num()].init(A);
    }
}

template <typename T, typename VectorType>
void spmv(const ThreadPartitionedSparseMatrix<T>& A, const VectorType& v, VectorType& res) {
    constexpr int64_t sizeofElem = (sizeof(T) + sizeof(A.matrices[0].innerIndexes[0]));
    RECORD_TIMER_PARAMS(A.nnz * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel
    {
        const typename ThreadPartitionedSparseMatrix<T>::SparseSubMatrix& localMat = A.matrices[omp_get_thread_num()];

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

template <typename T>
ThreadPartitionedSparseMatrix<T>::~ThreadPartitionedSparseMatrix() {
}

template <typename T>
struct ThreadPartitionedSparseMatrix<T>::SparseSubMatrix {
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

        rowStart = findPost(A, tid, numThreads);
        rowEnd = findPost(A, tid + 1, numThreads);

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

    template <typename other_t>
    static int findPost(const Eigen::SparseMatrix<other_t, MAJOR>& A, int blockId, int numBlocks) {
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

    int rowsGlobal;
    int rowStart;
    int rowEnd;
    int threadOwner;

    std::vector<int> innerIndexes;
    std::vector<int> outerIndexes;
    std::vector<T> data;
};

// void initFewSparseMatrices(const Eigen::SparseMatrix<double, MAJOR>& As,
//                            ThreadPartitionedSparseMatrix<double>& doubleMat,
//                            ThreadPartitionedSparseMatrix<float>& floatMat >) {
// }

template <typename T>
struct GreedyThreadPartitionedSparseMatrix<T>::GreedySparseSubMatrix {
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

        rowStart = findPost(A, tid, numThreads);
        rowEnd = findPost(A, tid + 1, numThreads);

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
                assert(pos != offsets.data() + offsetsCount);
                const int offsetIndex = (pos - offsets.data());
                assert(offsetIndex < offsetsCount);
                offsetInnerIndexes[j - dataOffset] = offsetIndex;
            }
        }
    }

    template <typename other_t>
    static int findPost(const Eigen::SparseMatrix<other_t, MAJOR>& A, int blockId, int numBlocks) {
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

    int rowsGlobal;
    int rowStart;
    int rowEnd;
    int threadOwner;

    std::vector<uint8_t> offsetInnerIndexes;
    std::vector<int> outerIndexes;
    std::vector<T> data;
};

static inline int mergeSorted_FIND_ANOTHER_PLACE_FOR_MY_IMPL_AND_DEF(const std::vector<std::array<int, 256>>& inArrays,
                                                                     std::vector<int> sizes,
                                                                     std::array<int, 256>& res) {
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

template <typename T>
template <typename other_t>
GreedyThreadPartitionedSparseMatrix<T>::GreedyThreadPartitionedSparseMatrix(
    const Eigen::SparseMatrix<other_t, MAJOR>& A)
    : nthr(omp_get_max_threads()), nnz(A.nonZeros()), matrices(nthr) {
    constexpr int64_t sizeofElem =
        (std::max(sizeof(T), sizeof(other_t)) + std::max(sizeof(int), sizeof(A.innerIndexPtr()[0])));
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

    const int offsetsCount =
        mergeSorted_FIND_ANOTHER_PLACE_FOR_MY_IMPL_AND_DEF(localOffsets, offsetsCountGlobal, colOffsets);

#pragma omp parallel
    {
        timer::commonTimer timerOmp("OMP section 2");
        matrices[omp_get_thread_num()].init(A, colOffsets, offsetsCount);
    }
}

template <typename T>
GreedyThreadPartitionedSparseMatrix<T>::~GreedyThreadPartitionedSparseMatrix() {
}

template <typename T, typename VectorType>
void spmv(const GreedyThreadPartitionedSparseMatrix<T>& A, const VectorType& v, VectorType& res) {
    constexpr int64_t sizeofElem = (sizeof(T) + sizeof(A.matrices[0].offsetInnerIndexes[0]));
    RECORD_TIMER_PARAMS(A.nnz * sizeofElem, timer::MeasureUnit::byte);
#pragma omp parallel
    {
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

template class GreedyThreadPartitionedSparseMatrix<double>;
template class GreedyThreadPartitionedSparseMatrix<float>;

template GreedyThreadPartitionedSparseMatrix<double>::GreedyThreadPartitionedSparseMatrix(
    const Eigen::SparseMatrix<double, MAJOR>& A);
template GreedyThreadPartitionedSparseMatrix<float>::GreedyThreadPartitionedSparseMatrix(
    const Eigen::SparseMatrix<double, MAJOR>& A);

template void spmv<double, Field3dBase<double>>(const GreedyThreadPartitionedSparseMatrix<double>& A,
                                                const Field3dBase<double>& v, Field3dBase<double>& res);
template void spmv<float, Field3dBase<float>>(const GreedyThreadPartitionedSparseMatrix<float>& A,
                                              const Field3dBase<float>& v, Field3dBase<float>& res);

template class ThreadPartitionedSparseMatrix<double>;
template class ThreadPartitionedSparseMatrix<float>;

template ThreadPartitionedSparseMatrix<double>::ThreadPartitionedSparseMatrix(
    const Eigen::SparseMatrix<double, MAJOR>& A);
template ThreadPartitionedSparseMatrix<float>::ThreadPartitionedSparseMatrix(
    const Eigen::SparseMatrix<double, MAJOR>& A);

template void spmv<double, Field3dBase<double>>(const ThreadPartitionedSparseMatrix<double>& A,
                                                const Field3dBase<double>& v, Field3dBase<double>& res);
template void spmv<float, Field3dBase<float>>(const ThreadPartitionedSparseMatrix<float>& A,
                                              const Field3dBase<float>& v, Field3dBase<float>& res);
