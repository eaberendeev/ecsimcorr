// Author: Evgeny Berendeev
// Email: evgeny.berendeev@gmail.com
// Copyright: (C) 2023, for licensing details see the LICENSE file

#include "util.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "sparse.h"
#include "timer.h"

// Функция, которая заполняет матрицу напрямую через внутренние массивы Eigen.
// Предполагается, что:
//   - trips отсортирован по (row, col)
//   - вектор trips не содержит дубликатов
//   - матрица хранится в формате RowMajor
void optimizedSetFromTriplets(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, const std::vector<Triplet> &trips) {
    RECORD_TIMER;

    int numRows = mat.rows();
    const int nnz = trips.size();

    // Устанавливаем размеры матрицы и выделяем ровно nnz элементов
    // mat.resize(numRows, numCols);
    mat.resizeNonZeros(nnz);

    // Получаем указатели на внутренние массивы
    int *outer = mat.outerIndexPtr();
    int *inner = mat.innerIndexPtr();
    double *values = mat.valuePtr();

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
#pragma omp parallel for schedule(guided)
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

Operator parallelSparseSum(const Operator &a, const Operator &b) {
    RECORD_TIMER;

    static_assert(a.IsRowMajor && b.IsRowMajor);
    assert(a.rows() == b.rows() && a.cols() == b.cols());
    assert(a.isCompressed() && b.isCompressed());

    const int rows = a.rows();

    // const Eigen::Map nnzA = a.innerNonZeros();
    const int *outerA = a.outerIndexPtr();
    const int *outerB = b.outerIndexPtr();

    const int *indA = a.innerIndexPtr();
    const int *indB = b.innerIndexPtr();

    std::vector<int> outerIndexes(rows + 1);
    outerIndexes[0] = 0.0;
    int nnz = 0;
    timer::commonTimer timerNNzCounter("nnz counter");
#pragma omp parallel for schedule(dynamic, 16 * 1024) reduction(+ : nnz)
    for (int i = 0; i < rows; ++i) {
        const int startA = outerA[i];
        const int endA = outerA[i + 1];
        const int startB = outerB[i];
        const int endB = outerB[i + 1];

        int nnzInRow = 0;

        int itA = startA;
        int itB = startB;
        while (itA != endA && itB != endB) {
            if (indA[itA] == indB[itB]) {
                itA += 1;
                itB += 1;
            } else if (indA[itA] < indB[itB]) {
                itA += 1;
            } else {
                itB += 1;
            }
            nnzInRow += 1;
        }

        nnzInRow += (endA - itA) + (endB - itB);
        outerIndexes[i + 1] = nnzInRow;
        nnz += nnzInRow;
    }
    timerNNzCounter.finish();

    for (int i = 1; i < rows + 1; ++i) {
        outerIndexes[i] += outerIndexes[i - 1];
    }

    Operator res(a.rows(), a.cols());
    res.resizeNonZeros(nnz);
    int *outerRes = res.outerIndexPtr();
    int *indRes = res.innerIndexPtr();
    // int* innerNnzRes = res.innerNonZeroPtr();

    const double *valuesA = a.valuePtr();
    const double *valuesB = b.valuePtr();
    double *valuesRes = res.valuePtr();

    outerRes[0] = 0;

    timer::commonTimer timerSummation("summation");
#pragma omp parallel for schedule(dynamic, 16 * 1024)
    for (int i = 0; i < rows; ++i) {
        outerRes[i + 1] = outerIndexes[i + 1];

        const int startA = outerA[i];
        const int endA = outerA[i + 1];
        const int startB = outerB[i];
        const int endB = outerB[i + 1];
        const int startRes = outerIndexes[i];

        int itA = startA;
        int itB = startB;
        int itRes = startRes;
        while (itA != endA && itB != endB) {
            if (indA[itA] == indB[itB]) {
                valuesRes[itRes] = valuesA[itA] + valuesB[itB];
                indRes[itRes] = indA[itA];
                itA += 1;
                itB += 1;
            } else if (indA[itA] < indB[itB]) {
                valuesRes[itRes] = valuesA[itA];
                indRes[itRes] = indA[itA];
                itA += 1;
            } else {
                valuesRes[itRes] = valuesB[itB];
                indRes[itRes] = indB[itB];
                itB += 1;
            }

            itRes += 1;
        }

        while (itA != endA) {
            valuesRes[itRes] = valuesA[itA];
            indRes[itRes] = indA[itA];
            itA += 1;
            itRes += 1;
        }
        while (itB != endB) {
            valuesRes[itRes] = valuesB[itB];
            indRes[itRes] = indB[itB];
            itB += 1;
            itRes += 1;
        }
    }
    timerSummation.finish();

    res.makeCompressed();
    return res;
}
