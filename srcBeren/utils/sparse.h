// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#pragma once

#ifndef SERVICE_H
#define SERVICE_H

#include <cassert>
#include <cmath>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>

#include "timer.h"
#include "types.h"
#include "util.h"

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
    ThreadPartitionedSparseMatrix(const Eigen::SparseMatrix<other_t, MAJOR>& A);

    ~ThreadPartitionedSparseMatrix();

    template <typename DataType, typename VectorType>
    friend void spmv(const ThreadPartitionedSparseMatrix<DataType>& A, const VectorType& v, VectorType& res);

    const int nthr;
    const int nnz;   // for timings only
    std::vector<SparseSubMatrix> matrices;
};

template <typename T>
struct ThreadPartitionedSparseMatrixView {
    ThreadPartitionedSparseMatrixView(const std::vector<int>& rowStartsIn, const std::vector<int>& rowEndsIn,
                                      const std::vector<std::vector<int>>& innerIndexesGlobIn,
                                      const std::vector<std::vector<int>>& outerIndexesGlobIn,
                                      const std::vector<std::vector<T>>& dataGlobIn)
        : rowStarts(rowStartsIn),
          rowEnds(rowEndsIn),
          innerIndexesGlob(innerIndexesGlobIn),
          outerIndexesGlob(outerIndexesGlobIn),
          dataGlob(dataGlobIn) {
    }

    template <typename VectorType>
    friend void spmv(const ThreadPartitionedSparseMatrixView& A, const VectorType& v, VectorType& res) {
        constexpr int64_t sizeofElem = (sizeof(T) + sizeof(A.innerIndexesGlob[0][0]));

#pragma omp parallel
        {
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

        const int maxThreads = omp_get_max_threads();
        rowStarts.resize(maxThreads);
        rowEnds.resize(maxThreads);
        innerIndexesGlob.resize(maxThreads);
        outerIndexesGlob.resize(maxThreads);
        dataGlob1.resize(maxThreads);
        dataGlob2.resize(maxThreads);

#pragma omp parallel num_threads(maxThreads)
        {
            int tid = omp_get_thread_num();
            std::vector<int>& innerIndexes = innerIndexesGlob[tid];
            std::vector<int>& outerIndexes = outerIndexesGlob[tid];
            std::vector<T1>& data1 = dataGlob1[tid];
            std::vector<T2>& data2 = dataGlob2[tid];

            const other_t* val = A.valuePtr();
            const int* inner = A.innerIndexPtr();
            const int* outer = A.outerIndexPtr();

            const int rowStart = findBestPos(A, tid, maxThreads);
            const int rowEnd = findBestPos(A, tid + 1, maxThreads);
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
        if constexpr (std::is_same_v<T, T1>) {
            return ThreadPartitionedSparseMatrixView<T>(rowStarts, rowEnds, innerIndexesGlob, outerIndexesGlob,
                                                        dataGlob1);
        }
        if constexpr (std::is_same_v<T, T2>)
            return ThreadPartitionedSparseMatrixView<T>(rowStarts, rowEnds, innerIndexesGlob, outerIndexesGlob,
                                                        dataGlob2);
        static_assert(std::is_same_v<T, T1> || std::is_same_v<T, T2>);
    }

    template <typename other_t>
    static int findBestPos(const Eigen::SparseMatrix<other_t, MAJOR>& A, int blockId, int numBlocks) {
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

    std::vector<int> rowStarts;
    std::vector<int> rowEnds;
    std::vector<std::vector<int>> innerIndexesGlob;
    std::vector<std::vector<int>> outerIndexesGlob;
    std::vector<std::vector<T1>> dataGlob1;
    std::vector<std::vector<T2>> dataGlob2;
};

template <typename T>
struct GreedyThreadPartitionedSparseMatrix {
   private:
    struct GreedySparseSubMatrix;

   public:
    template <typename other_t>
    GreedyThreadPartitionedSparseMatrix(const Eigen::SparseMatrix<other_t, MAJOR>& A);

    ~GreedyThreadPartitionedSparseMatrix();

    template <typename DataType, typename VectorType>
    friend void spmv(const GreedyThreadPartitionedSparseMatrix<DataType>& A, const VectorType& v, VectorType& res);

    std::array<int, 256> colOffsets;
    int offsetsCount;
    const int nthr;
    const int nnz;   // for timings only
    std::vector<GreedySparseSubMatrix> matrices;
};

template <typename VectorType>
inline void spmv(const Operator& A, const VectorType& v, VectorType& res) {
    RECORD_TIMER_PARAMS(A.nonZeros() * (sizeof(double) + sizeof(A.innerIndexPtr()[0])), timer::MeasureUnit::byte);
    int rows = A.rows();

    const double* val = A.valuePtr();
    const int* inner = A.innerIndexPtr();
    const int* outer = A.outerIndexPtr();

    std::vector<int> diffs;
    std::map<int, int> colsCount;
    std::map<int, int> singleElemsOffsets;
    std::map<std::vector<int>, int> colsOffsets;

#pragma omp parallel
    {
        timer::commonTimer timerOmp("OMP section");
#pragma omp for schedule(dynamic, 16 * 1024)
        for (int i = 0; i < rows; ++i) {
            double sum = 0;
#pragma omp simd
            for (int j = outer[i]; j < outer[i + 1]; ++j) {
                sum += val[j] * v[inner[j]];
            }
            res[i] = sum;
        }
    }

    int totalSize = 0;
    for (const auto& it : colsCount) {
        totalSize += it.second * it.first;
    }
}

// Структура для хранения элемента разреженной матрицы
class Triplet {
   public:
    Triplet(int r, int c, double v) : _row(r), _col(c), _value(v) {
    }
    Triplet() : _row(0), _col(0), _value(0) {
    }
    const int& row() const noexcept {
        return _row;
    }
    const int& col() const noexcept {
        return _col;
    }
    const double& value() const noexcept {
        return _value;
    }
    int& row() {
        return _row;
    }
    int& col() {
        return _col;
    }
    double& value() {
        return _value;
    }
    // Triplet(Triplet&& other) noexcept = default;
    // Triplet& operator=(Triplet&& other) noexcept = default;
    bool operator<(const Triplet& other) const {
        return std::tie(_row, _col) < std::tie(other.row(), other.col());
    }

   private:
    int _row;
    int _col;
    double _value;
};

inline bool compareTriplets(const Triplet& a, const Triplet& b) {
    return std::tie(a.row(), a.col()) < std::tie(b.row(), b.col());
}

void optimizedSetFromSortedTriplets(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                                    const std::vector<Triplet>& trips);

inline Operator parallel_sparse_addition(const Operator& A, double alpha, const Operator& B, double beta) {
    int rows = A.rows();
    int cols = A.cols();
    int num_threads = omp_get_max_threads();

    // Проверка на совместимость матриц
    if (B.rows() != rows || B.cols() != cols) {
        throw std::invalid_argument("Размеры матриц A и B должны совпадать");
    }

    // Векторы для хранения данных каждого потока
    std::vector<std::vector<double>> thread_values(num_threads);
    std::vector<std::vector<int>> thread_col_indices(num_threads);
    std::vector<std::vector<int>> thread_row_pointers(num_threads, std::vector<int>(rows + 1, 0));

// Параллельная обработка
#pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();

        std::vector<double> local_values;
        std::vector<int> local_col_indices;
        std::vector<int> local_row_pointers(rows + 1, 0);

        // Резервирование памяти
        local_values.reserve(A.nonZeros() / num_threads);
        local_col_indices.reserve(A.nonZeros() / num_threads);

#pragma omp for schedule(dynamic, 512)
        for (int i = 0; i < rows; ++i) {
            Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_a(A, i);
            Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_b(B, i);

            int row_start = local_values.size();
            while (it_a || it_b) {
                int col;
                double value = 0.0;

                if (it_a && (!it_b || it_a.col() < it_b.col())) {
                    col = it_a.col();
                    value = alpha * it_a.value();
                    ++it_a;
                } else if (it_b && (!it_a || it_b.col() < it_a.col())) {
                    col = it_b.col();
                    value = beta * it_b.value();
                    ++it_b;
                } else {
                    col = it_a.col();
                    value = alpha * it_a.value() + beta * it_b.value();
                    ++it_a;
                    ++it_b;
                }

                if (value != 0.0) {
                    local_values.push_back(value);
                    local_col_indices.push_back(col);
                }
            }
            local_row_pointers[i + 1] = local_values.size() - row_start;
        }

        // Преобразование локальных указателей строк в абсолютные индексы
        for (int i = 1; i <= rows; ++i) {
            local_row_pointers[i] += local_row_pointers[i - 1];
        }

        // Сохраняем результаты потока
        thread_values[thread_id] = std::move(local_values);
        thread_col_indices[thread_id] = std::move(local_col_indices);
        thread_row_pointers[thread_id] = std::move(local_row_pointers);
    }

    // Считаем общее количество ненулевых элементов
    int total_nnz = 0;
    for (const auto& vals : thread_values) {
        total_nnz += vals.size();
    }
    Operator C(rows, cols);
    // Резервируем память для итоговой матрицы
    C.reserve(total_nnz);

    double* C_values = C.valuePtr();
    int* C_inner_indices = C.innerIndexPtr();
    int* C_outer_indices = C.outerIndexPtr();

    // Слияние данных всех потоков
    std::vector<int> row_start(rows + 1, 0);   // Указатели начала строк в глобальной матрице
    int nnz_offset = 0;
    for (int i = 0; i < rows; ++i) {
        for (int t = 0; t < num_threads; ++t) {
            if (thread_row_pointers[t][i + 1] > thread_row_pointers[t][i]) {
                row_start[i] = nnz_offset;
                break;
            }
        }
        for (int t = 0; t < num_threads; ++t) {
            int local_start = thread_row_pointers[t][i];
            int local_end = thread_row_pointers[t][i + 1];
            if (local_end > local_start) {
                int chunk_size = local_end - local_start;
                std::copy(thread_values[t].begin() + local_start, thread_values[t].begin() + local_end,
                          C_values + nnz_offset);
                std::copy(thread_col_indices[t].begin() + local_start, thread_col_indices[t].begin() + local_end,
                          C_inner_indices + nnz_offset);
                nnz_offset += chunk_size;
            }
        }
        C_outer_indices[i + 1] = nnz_offset;
    }

    // Финализируем матрицу
    C_outer_indices[rows] = nnz_offset;
    C.finalize();
    return C;
}

inline Operator parallel_sparse_addition2(const Operator& A, double alpha, const Operator& B, double beta) {
    int rows = A.rows();
    int cols = A.cols();
    int num_threads = omp_get_max_threads();

    // Проверка совместимости размеров
    if (B.rows() != rows || B.cols() != cols) {
        throw std::invalid_argument("Размеры матриц A и B должны совпадать");
    }

    // Векторы для хранения данных каждого потока
    std::vector<std::vector<double>> thread_values(num_threads);
    std::vector<std::vector<int>> thread_col_indices(num_threads);
    // Каждый поток хранит вектор размером rows+1, где для строки i:
    // thread_row_pointers[t][i+1] - thread_row_pointers[t][i] = число ненулевых
    // элементов в строке i, обработанных этим потоком.
    std::vector<std::vector<int>> thread_row_pointers(num_threads, std::vector<int>(rows + 1, 0));

    // Параллельная обработка строк
#pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();

        std::vector<double> local_values;
        std::vector<int> local_col_indices;
        std::vector<int> local_row_pointers(rows + 1, 0);

        // Резервирование памяти (может быть уточнено, если известна статистика)
        local_values.reserve(A.nonZeros() / num_threads);
        local_col_indices.reserve(A.nonZeros() / num_threads);

#pragma omp for schedule(dynamic, 512)
        for (int i = 0; i < rows; ++i) {
            Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_a(A, i);
            Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_b(B, i);

            int row_start = local_values.size();
            // Слияние строк из A и B (предполагается, что строки отсортированы
            // по столбцам)
            while (it_a || it_b) {
                int col;
                double value = 0.0;

                if (it_a && (!it_b || it_a.col() < it_b.col())) {
                    col = it_a.col();
                    value = alpha * it_a.value();
                    ++it_a;
                } else if (it_b && (!it_a || it_b.col() < it_a.col())) {
                    col = it_b.col();
                    value = beta * it_b.value();
                    ++it_b;
                } else {
                    col = it_a.col();
                    value = alpha * it_a.value() + beta * it_b.value();
                    ++it_a;
                    ++it_b;
                }

                if (value != 0.0) {
                    local_values.push_back(value);
                    local_col_indices.push_back(col);
                }
            }
            // Запоминаем количество внесённых элементов для строки i
            local_row_pointers[i + 1] = local_values.size() - row_start;
        }

        // Преобразование локальных указателей строк в абсолютные индексы
        for (int i = 1; i <= rows; ++i) {
            local_row_pointers[i] += local_row_pointers[i - 1];
        }

        // Сохраняем результаты потока
        thread_values[thread_id] = std::move(local_values);
        thread_col_indices[thread_id] = std::move(local_col_indices);
        thread_row_pointers[thread_id] = std::move(local_row_pointers);
    }

    // Вычисляем для каждой строки общее число ненулевых элементов,
    // суммируя вклады от всех потоков.
    std::vector<int> global_row_counts(rows, 0);
    for (int i = 0; i < rows; ++i) {
        for (int t = 0; t < num_threads; ++t) {
            global_row_counts[i] += thread_row_pointers[t][i + 1] - thread_row_pointers[t][i];
        }
    }

    // Вычисляем глобальные указатели строк (префиксное суммирование)
    std::vector<int> C_outer_indices(rows + 1, 0);
    for (int i = 0; i < rows; ++i) {
        C_outer_indices[i + 1] = C_outer_indices[i] + global_row_counts[i];
    }
    int total_nnz = C_outer_indices[rows];

    // Инициализируем итоговую матрицу и резервируем память
    Operator C(rows, cols);
    C.reserve(total_nnz);
    double* C_values = C.valuePtr();
    int* C_inner_indices = C.innerIndexPtr();
    int* C_outer_ptr = C.outerIndexPtr();

    // Копируем рассчитанные указатели строк в итоговую матрицу
    std::copy(C_outer_indices.begin(), C_outer_indices.end(), C_outer_ptr);

    // Параллельное слияние данных по строкам
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < rows; ++i) {
        int offset = C_outer_indices[i];
        // Для каждой строки копируем данные из локальных векторов каждого
        // потока
        for (int t = 0; t < num_threads; ++t) {
            int local_start = thread_row_pointers[t][i];
            int local_end = thread_row_pointers[t][i + 1];
            int count = local_end - local_start;
            if (count > 0) {
                std::copy(thread_values[t].begin() + local_start, thread_values[t].begin() + local_end,
                          C_values + offset);
                std::copy(thread_col_indices[t].begin() + local_start, thread_col_indices[t].begin() + local_end,
                          C_inner_indices + offset);
                offset += count;
            }
        }
    }

    C.finalize();
    return C;
}

inline Operator parallel_sparse_addition3(const Operator& A, double alpha, const Operator& B, double beta) {
    int rows = A.rows();
    int cols = A.cols();

    if (B.rows() != rows || B.cols() != cols) {
        throw std::invalid_argument("Matrix dimensions must match");
    }

    // Первый проход: вычисление количества ненулевых элементов для каждой
    // строки
    std::vector<int> nnz_per_row(rows, 0);

#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < rows; ++i) {
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_a(A, i);
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_b(B, i);
        int count = 0;

        while (it_a || it_b) {
            int col_a = (it_a) ? it_a.col() : cols;
            int col_b = (it_b) ? it_b.col() : cols;

            if (col_a < col_b) {
                double val = alpha * it_a.value();
                if (val != 0.0)
                    ++count;   // Проверяем val, а не alpha != 0
                ++it_a;
            } else if (col_b < col_a) {
                double val = beta * it_b.value();
                if (val != 0.0)
                    ++count;   // Проверяем val, а не beta != 0
                ++it_b;
            } else {
                double val = alpha * it_a.value() + beta * it_b.value();
                if (val != 0.0)
                    ++count;
                ++it_a;
                ++it_b;
            }
        }

        nnz_per_row[i] = count;
    }

    // Инициализация матрицы C с резервированием памяти
    Operator C(rows, cols);
    C.reserve(nnz_per_row);

    double* C_values = C.valuePtr();
    int* C_inner = C.innerIndexPtr();
    int* C_outer = C.outerIndexPtr();

    // Заполнение C_outer на основе nnz_per_row
    C_outer[0] = 0;
    for (int i = 0; i < rows; ++i) {
        C_outer[i + 1] = C_outer[i] + nnz_per_row[i];
    }

// Второй проход: заполнение данных матрицы C
#pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < rows; ++i) {
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_a(A, i);
        Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_b(B, i);
        int pos = C_outer[i];

        while (it_a || it_b) {
            int col_a = (it_a) ? it_a.col() : std::numeric_limits<int>::max();
            int col_b = (it_b) ? it_b.col() : std::numeric_limits<int>::max();

            if (col_a < col_b) {
                double val = alpha * it_a.value();
                if (val != 0.0) {
                    C_values[pos] = val;
                    C_inner[pos] = col_a;
                    pos++;
                }
                if (it_a)
                    ++it_a;
            } else if (col_b < col_a) {
                double val = beta * it_b.value();
                if (val != 0.0) {
                    C_values[pos] = val;
                    C_inner[pos] = col_b;
                    pos++;
                }
                if (it_b)
                    ++it_b;
            } else {
                double val = alpha * it_a.value() + beta * it_b.value();
                if (val != 0.0) {
                    C_values[pos] = val;
                    C_inner[pos] = col_a;
                    pos++;
                }
                if (it_a)
                    ++it_a;
                if (it_b)
                    ++it_b;
            }
        }
    }

    C.finalize();
    return C;
}

Operator parallelSparseSum(const Operator& a, const Operator& b);

bool checkMatrixPortraitCoincidence(const Operator& a, const Operator& b);

// for debug purposes only
void checkMatrixCoincidence(const Operator& a, const Operator& b, const double relTolerance);

#endif   // SERVICE_H
