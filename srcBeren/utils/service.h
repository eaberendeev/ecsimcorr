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

#include "util.h"

// general indexing routine (row major)
// inline constexpr int ind(int x, int y, int z, int c, int Nx, int Ny, int Nz,
//                          int Nc) {
//     return (c + Nc * (z + Nz * (y + Ny * x)) + 0 * Nx);
// }

bool create_directory(const std::string& path);
std::vector<std::string> split_string(const std::string& s, const char delim);

inline int int_value(const double dblValue) {
    return static_cast<int>(dblValue < 0 ? dblValue - 0.5 : dblValue + 0.5);
}

template <typename T>
inline int sign(const T val) {
    return (T(0) < val) - (val < T(0));
}

inline std::string to_string(const int i, const int len) {
    std::stringstream ss;
    ss << std::setw(len) << std::setfill('0') << i;
    std::string s = ss.str();
    return s;
}

inline Vector3i cast_to_int(const Vector3d& value) {
    return Vector3i(static_cast<int>(value(0)), static_cast<int>(value(1)),
                    static_cast<int>(value(2)));
}

inline Vector3i vector_to_vector3i(const std::vector<int>& vec) {
    assert(vec.size() == 3 &&
           "Std::vector has invalid size. Convertation to Vecto3d: FAILED");
    return Vector3i(vec[0], vec[1], vec[2]);
}
inline Vector3d vector_to_vector3d(const std::vector<double>& vec) {
    assert(vec.size() == 3 &&
           "Std::vector has invalid size. Convertation to Vecto3d: FAILED");
    return Vector3d(vec[0], vec[1], vec[2]);
}

inline std::vector<double> string_to_doubles(std::vector<std::string> words) {
    std::vector<double> doubles;
    for (const auto& word : words) {
        doubles.push_back(std::stod(word));
    }
    return doubles;
}

template <typename Func>
std::pair<double, double> find_maximum_universal(Func f, double x_min,
                                                 double x_max, int intervals) {
    // if (intervals <= 0 || x_min >= x_max) {
    //     throw std::invalid_argument("Invalid input parameters");
    // }

    double step = (x_max - x_min) / intervals;
    double max_value = f(x_min);
    double x_optimal = x_min;

    for (int i = 1; i <= intervals; ++i) {
        double x = x_min + i * step;
        double current_value = f(x);
        if (current_value > max_value) {
            max_value = current_value;
            x_optimal = x;
        }
    }

    return {max_value, x_optimal};
}

inline int hash(const std::string& key, int tableSize) {
    int hashVal = 0;

    for (size_t i = 0; i < key.length(); i++) {
        hashVal = 37 * hashVal + key[i];
    }

    hashVal %= tableSize;

    if (hashVal < 0)
        hashVal += tableSize;

    return hashVal;
}

struct pair_hash {
    size_t operator()(const std::pair<int, int>& p) const {
        return static_cast<size_t>(p.first) << 32 | p.second;
    }
};

// Структура для хранения элемента разреженной матрицы
class Triplet {

public:
    Triplet(int r, int c, double v) : _row(r), _col(c), _value(v) {}
    Triplet() : _row(0), _col(0), _value(0) {}
    const int& row() const noexcept { return _row; }
    const int& col() const noexcept { return _col; }
    const double& value() const noexcept { return _value; }
    int& row() {return _row;}
    int& col() {return _col;}
    double& value() {return _value;}
    // Triplet(Triplet&& other) noexcept = default;
    //Triplet& operator=(Triplet&& other) noexcept = default;
    bool operator < (const Triplet& other) const {
        return std::tie(_row, _col) < std::tie(other.row(), other.col());
    }
    private:
     int _row;
     int _col;
     double _value;
};

void optimizedSetFromTriplets(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat,
                              const std::vector<Triplet>& trips);

template <typename T>
T get_checked(const nlohmann::json& j, const std::string& key) {
    if (!j.contains(key)) {
        throw std::runtime_error("Missing key: " + key);
    }
    return j[key].get<T>();
}

// Специализации для более информативных сообщений
template <>
inline std::string get_checked<std::string>(const nlohmann::json& j,
                                     const std::string& key) {
    if (!j.contains(key)) {
        throw std::runtime_error("Missing key: " + key);
    }

    if (!j[key].is_string()) {
        throw std::runtime_error("Key '" + key + "' is not a string");
    }

    return j[key].get<std::string>();
}

template <>
inline int get_checked<int>(const nlohmann::json& j, const std::string& key) {
    if (!j.contains(key)) {
        throw std::runtime_error("Missing key: " + key);
    }

    if (!j[key].is_number_integer()) {
        throw std::runtime_error("Key '" + key + "' is not an integer");
    }

    return j[key].get<int>();
}


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

#endif   // SERVICE_H
