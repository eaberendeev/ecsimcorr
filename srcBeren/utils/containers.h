#ifndef CONTAINERS_H_
#define CONTAINERS_H_
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Vec.h"
#include "util.h"
//#define NDEBUG 0
// Add or subtract vector b to vector a element-wise
// Assumes a and b are the same size
// Enable only for arithmetic types
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());

    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<T>());

    return a;
}
template <typename T, std::enable_if_t<std::is_arithmetic<T>::value>* = nullptr>
std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size() == b.size());

    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<T>());

    return a;
}

template <typename T>
class Array3D {
   public:
    Array3D(int n1, int n2, int n3) { resize(n1, n2, n3); }
    Array3D(const int3& nn) { resize(nn(0), nn(1), nn(2)); }
    Array3D(Array3D&& other) : _data{other._data}, _size{other._size} {}

    Array3D() { _size = int3(0, 0, 0); }

    void resize(int n1, int n2, int n3) {
        _data.resize(n1 * n2 * n3);
        _size = int3(n1, n2, n3);
    }
    void resize(const int3& nn) {
        _data.resize(nn.total_size());
        _size = nn;
    }
    void setZero() {
        for (auto i = 0; i < capacity(); ++i) {
            _data[i] = 0.;
        }
    }

    Array3D<T>& operator=(Array3D<T>&& array3D) {
        // Self-assignment check
        if (this == &array3D)
            return *this;
        // Checking the conformity of dimensions
        assert(_size == array3D._size);
        _data = array3D._data;
        _size = array3D._size;
        return (*this);
    }

    bool checkIndex(int i, int j, int k) const {
        assert(i < _size(0) && j < _size(1) && k < _size(2));
        assert(i >= 0 && j >= 0 && k >= 0);
        return true;
    }

    T& operator()(int i, int j, int k) {
        assert(checkIndex(i, j, k));
        return _data[i * _size(1) * _size(2) + j * _size(2) + k];
    }

    const T& operator()(int i, int j, int k) const {
        assert(checkIndex(i, j, k));
        return _data[i * _size(1) * _size(2) + j * _size(2) + k];
    }
    T& operator()(int i) { return _data[i]; }
    const T& operator()(int i) const { return _data[i]; }
    std::vector<T>& data() { return _data; }
    const std::vector<T>& data() const { return _data; }

    int3 size() const { return _size; }
    int capacity() const { return _size.total_size(); }

   private:
    std::vector<T> _data;
    int3 _size;
};

template <typename T>
class Array2D {
   public:
    Array2D(int n1, int n2) { resize(n1, n2); }
    Array2D(const int2& nn) { resize(nn(0), nn(1)); }
    Array2D(Array2D&& other) : _data{other._data}, _size{other._size} {}
    Array2D() { _size = int2(0, 0, 0); }

    void resize(int n1, int n2) {
        _data.resize(n1 * n2);
        _size = int2(n1, n2);
    }
    void resize(const int2& nn) {
        _data.resize(nn.total_size());
        _size = nn;
    }
    void setZero() {
        for (auto i = 0; i < capacity(); ++i) {
            _data[i] = 0.;
        }
    }
    Array2D<T>& operator=(Array2D<T>&& other) {
        // Self-assignment check
        if (this == &other)
            return *this;
        // Checking the conformity of dimensions
        assert(_size == other._size);
        _data = other._data;
        _size = other._size;
        return (*this);
    }
    bool checkIndex(int i, int j) const {
        assert(i < _size(0) && j < _size(1));
        assert(i >= 0 && j >= 0);
        return true;
    }
    T& operator()(int i, int j) {
        assert(checkIndex(i, j));
        return _data[i * _size(1) + j];
    }
    const T& operator()(int i, int j) const {
        assert(checkIndex(i, j));
        return _data[i * _size(1) + j];
    }
    T& operator()(int i) { return _data[i]; }
    const T& operator()(int i) const { return _data[i]; }
    std::vector<T>& data() { return _data; }
    const std::vector<T>& data() const { return _data; }

    int2 size() const { return _size; }
    int capacity() const { return _size.total_size(); }

   private:
    std::vector<T> _data;
    int2 _size;
};

class Field3d {
   public:
    Field3d(const int& n) { resize(n); }

    Field3d(int n1, int n2, int n3, int d) { resize(n1, n2, n3, d); }
    Field3d(const int3& nn, int d) { resize(nn, d); }
    Field3d(const Field3d& other) : _size{other._size}, _nd{other._nd} {
        resize(_size, _nd);
#pragma omp parallel for simd
        for (long i = 0; i < _data.size(); ++i) {
            _data[i] = other._data[i];
        }
    }
    // move constructor
    Field3d(Field3d&& other) noexcept
        : _data(std::move(other._data)), _size(other._size), _nd(other._nd) {
        other._size = int3(0, 0, 0);
        other._nd = 0;
    }
    Field3d() { _size = int3(0, 0, 0); }
    // Add static Zero constructor
    static Field3d Zero(const int3& size, int nd) {
        Field3d result(size, nd);
#pragma omp parallel for simd
        for (int i = 0; i < result.capacity(); ++i) {
            result._data[i] = 0.0;
        }
        return result;
    }
    static Field3d Zero(int n) {
        Field3d result(n);
#pragma omp parallel for simd
        for (long i = 0; i < result.capacity(); ++i) {
            result._data[i] = 0.0;
        }
        return result;
    }
    void resize(int n1, int n2, int n3, int d) {
        _data.resize(n1 * n2 * n3 * d);
        _size = int3(n1, n2, n3);
        _nd = d;
    }
    void resize(const int3& nn, int d) {
        _data.resize(nn.total_size() * d);
        _size = nn;
        _nd = d;
    }
    void resize(int n){
        _data.resize(n);
        _size = int3(n, 1, 1);
        _nd = 1;
    }
    void setZero() {
#pragma omp parallel for simd
        for (auto i = 0; i < capacity(); ++i) {
            _data[i] = 0.;
        }
    }
    bool checkIndex(int i, int j, int k, int d) const {
        bool true1 = i < _size(0) && j < _size(1) && k < _size(2) && d < _nd;
        bool true2 = i >= 0 && j >= 0 && k >= 0 && d >= 0;
        return true1 && true2;
    }
    Field3d& operator=(const Field3d& other) {
        // Self-assignment check
        if (this == &other)
            return *this;
        // Checking the conformity of dimensions
        assert(_size == other._size && _nd == other._nd);

        //_data = other._data;
#pragma omp parallel for simd
        for(long i = 0; i < _data.size(); ++i) {
            _data[i] = other._data[i];
        }
        return (*this);
    }
    // move assignment
    Field3d& operator=(Field3d&& other) noexcept {
        if (this != &other) {
            _data = std::move(other._data);
            _size = other._size;
            _nd = other._nd;
            other._size = int3(0, 0, 0);
            other._nd = 0;
        }
        return *this;
    }
    double& operator()(int i, int j, int k, int d) {
        assert(checkIndex(i, j, k, d));
        return _data[d + _nd * (i * _size(1) * _size(2) + j * _size(2) + k)];
    }
    const double& operator()(int i, int j, int k, int d) const {
        assert(checkIndex(i, j, k, d));
        return _data[d + _nd * (i * _size(1) * _size(2) + j * _size(2) + k)];
    }
    double& operator()(int i) { return _data[i]; }
    const double& operator()(int i) const { return _data[i]; }
    double& operator[](int i) { return _data[i]; }
    const double& operator[](int i) const { return _data[i]; }
    int nd() const { return _nd; }
    int3 sizes() const { return _size; }
    int capacity() const { return _nd * _size.total_size(); }
    size_t size() const { return (size_t)_nd * _size.total_size(); }
    Eigen::VectorXd& data() { return _data; }
    const Eigen::VectorXd& data() const { return _data; }

    // Scalar multiplication
    Field3d& operator*=(const double alpha) {
#pragma omp parallel for simd
        for (int i = 0; i < capacity(); ++i) {
            _data[i] *= alpha;
        }
        return *this;
    }

    // Addition and subtraction
    Field3d& operator+=(const Field3d& other) {
        assert(capacity() == other.capacity());
#pragma omp parallel for simd
        for (int i = 0; i < capacity(); ++i) {
            _data[i] += other._data[i];
        }
        return *this;
    }

    Field3d& operator-=(const Field3d& other) {
        assert(capacity() == other.capacity());
#pragma omp parallel for simd
        for (int i = 0; i < capacity(); ++i) {
            _data[i] -= other._data[i];
        }
        return *this;
    }


    // dot product
    double dot(const Field3d& other) const {
        assert(capacity() == other.capacity());
        double result = 0.0;
#pragma omp parallel for simd reduction(+ : result)
        for (int i = 0; i < capacity(); ++i) {
            result += _data[i] * other._data[i];
        }
        return result;
    }

    // squared norm
    double squaredNorm() const {
        double result = 0.0;
#pragma omp parallel for simd reduction(+ : result)
     for (int i = 0; i < capacity(); ++i) {
        result += _data[i] * _data[i];
    }
    return result;
    }

    double norm() const { return std::sqrt(squaredNorm()); }

    friend Field3d operator*(const Field3d& field, const double alpha) {
        Field3d result(field);
        result *= alpha;
        return result;
    }
    friend Field3d operator*(const double alpha, const Field3d& field) {
        return field * alpha;
    }
    friend Field3d operator*(const Operator& A, const Field3d& field) {
        int rows = A.rows();
        Field3d res = Field3d(field);
#pragma omp parallel for simd schedule(dynamic, 256)
        for (int row = 0; row < rows; row++) {
            double res_row = 0.0;
            for (Operator::InnerIterator it(A, row); it; ++it) {
                res_row += it.value() * field(it.col());
            }
            res(row) = res_row;
        }

        return res;
    }
    friend Field3d operator+(const Field3d& a, const Field3d& b) {
        Field3d result(a);
        result += b;
        return result;
    }

    friend Field3d operator-(const Field3d& a, const Field3d& b) {
        Field3d result(a);
        result -= b;
        return result;
    }

    // Add iterator support
    using iterator = typename Eigen::VectorXd::iterator;
    using const_iterator = typename Eigen::VectorXd::const_iterator;

    iterator begin() { return _data.begin(); }
    iterator end() { return _data.end(); }
    const_iterator begin() const { return _data.begin(); }
    const_iterator end() const { return _data.end(); }
    const_iterator cbegin() const { return _data.cbegin(); }
    const_iterator cend() const { return _data.cend(); }

   private:
    Eigen::VectorXd _data;
    int3 _size;
    int _nd;
};

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
    std::vector<std::vector<int>> thread_row_pointers(
        num_threads, std::vector<int>(rows + 1, 0));

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
                    value = alpha*it_a.value();
                    ++it_a;
                } else if (it_b && (!it_a || it_b.col() < it_a.col())) {
                    col = it_b.col();
                    value = beta*it_b.value();
                    ++it_b;
                } else {
                    col = it_a.col();
                    value = alpha*it_a.value() + beta*it_b.value();
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
    std::vector<int> row_start(
        rows + 1, 0);   // Указатели начала строк в глобальной матрице
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
                std::copy(thread_values[t].begin() + local_start,
                          thread_values[t].begin() + local_end,
                          C_values + nnz_offset);
                std::copy(thread_col_indices[t].begin() + local_start,
                          thread_col_indices[t].begin() + local_end,
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
    // thread_row_pointers[t][i+1] - thread_row_pointers[t][i] = число ненулевых элементов в строке i, обработанных этим потоком.
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
            // Слияние строк из A и B (предполагается, что строки отсортированы по столбцам)
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
        // Для каждой строки копируем данные из локальных векторов каждого потока
        for (int t = 0; t < num_threads; ++t) {
            int local_start = thread_row_pointers[t][i];
            int local_end = thread_row_pointers[t][i + 1];
            int count = local_end - local_start;
            if (count > 0) {
                std::copy(thread_values[t].begin() + local_start,
                          thread_values[t].begin() + local_end,
                          C_values + offset);
                std::copy(thread_col_indices[t].begin() + local_start,
                          thread_col_indices[t].begin() + local_end,
                          C_inner_indices + offset);
                offset += count;
            }
        }
    }

    C.finalize();
    return C;
}


inline Operator parallel_sparse_addition3(const Operator& A, double alpha,
                                  const Operator& B, double beta) {
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

#endif
