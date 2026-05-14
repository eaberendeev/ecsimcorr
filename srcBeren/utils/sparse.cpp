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

// Функция, которая заполняет матрицу напрямую через внутренние массивы Eigen.
// Предполагается, что:
//   - trips отсортирован по (row, col)
//   - вектор trips не содержит дубликатов
//   - матрица хранится в формате RowMajor
void optimizedSetFromTriplets(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat, const std::vector<Triplet>& trips) {
    int numRows = mat.rows();
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
