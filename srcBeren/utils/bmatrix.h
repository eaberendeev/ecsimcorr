#pragma once

#include <assert.h>

#include <unordered_map>
#include <vector>

#include "util.h"

typedef std::unordered_map<int, double> IndexMap;

class BMatrix {
   public:
    BMatrix(size_t size) : data(size) {}
    BMatrix() {}
    BMatrix(const BMatrix& other) : data(other.data) {}
    BMatrix(BMatrix&& other) : data(other.data) {}
    BMatrix(const Operator& mat){
        data.resize(mat.rows());
        for (int k = 0; k < mat.outerSize(); ++k) {
            for (Operator::InnerIterator it(mat, k); it; ++it) {
                data[it.row()][it.col()] = it.value();
            }
        }
    }
    // Move assignment operator
    BMatrix& operator=(BMatrix&& other) noexcept {
        if (this != &other) {
            data = std::move(other.data);
        }
        return *this;
    }
    void reserve(size_t size) { data.reserve(size); }
    void resize(size_t size) { data.resize(size); }
    // Assignment operator
    BMatrix& operator=(const BMatrix& other) {
        if (this != &other) {
            data = other.data;
        }
        return *this;
    }

    // Addition assignment operator
    BMatrix& operator+=(const BMatrix& other) {
        for (size_t i = 0; i < data.size(); i++) {
            for (const auto& [j, val] : other.data[i]) {
                data[i][j] += val;
            }
        }
        return *this;
    }

    // Subtraction assignment operator
    BMatrix& operator-=(const BMatrix& other) {
        for (size_t i = 0; i < data.size(); i++) {
            for (const auto& [j, val] : other.data[i]) {
                data[i][j] -= val;
            }
        }
        return *this;
    }

    // Matrix multiplication
    BMatrix operator*(const BMatrix& other) const {
        BMatrix result(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            for (const auto& [k, v] : data[i]) {
                for (const auto& [j, val] : other.data[k]) {
                    result.data[i][j] += v * val;
                }
            }
        }
        return result;
    }

    // Matrix addition
    BMatrix operator+(const BMatrix& other) const {
        BMatrix result(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            result.data[i] = data[i];
            for (const auto& [j, val] : other.data[i]) {
                result.data[i][j] += val;
            }
        }
        return result;
    }

    // Matrix subtraction
    BMatrix operator-(const BMatrix& other) const {
        BMatrix result(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            result.data[i] = data[i];
            for (const auto& [j, val] : other.data[i]) {
                result.data[i][j] -= val;
            }
        }
        return result;
    }

    // Matrix-vector multiplication
    std::vector<double> operator*(const std::vector<double>& vec) const {
        std::vector<double> result(data.size(), 0.0);
        for (size_t i = 0; i < data.size(); i++) {
            for (const auto& [j, val] : data[i]) {
                result[i] += val * vec[j];
            }
        }
        return result;
    }

    // Matrix-Field multiplication
    Field operator*(const Field& vec) const {
        Field result = Field::Zero(data.size());
        for (size_t i = 0; i < data.size(); i++) {
            for (const auto& [j, val] : data[i]) {
                result[i] += val * vec[j];
            }
        }
        return result;
    }
    // Scalar multiplication (scalar * matrix)
    friend BMatrix operator*(double scalar, const BMatrix& mat) {
        BMatrix result(mat.size());
        for (size_t i = 0; i < mat.size(); i++) {
            for (const auto& [j, val] : mat.data[i]) {
                result.data[i][j] = scalar * val;
            }
        }
        return result;
    }
    // Scalar multiplication (matrix * scalar)
    BMatrix operator*(double scalar) { return scalar * (*this); }
    // Scalar multiplication assignment
    BMatrix& operator*=(double scalar) {
        for (auto& row : data) {
            for (auto& [_, val] : row) {
                val *= scalar;
            }
        }
        return *this;
    }
    // Operator[] for accessing rows
    IndexMap& operator[](size_t index) {
        assert(index < data.size());
        return data[index];
    }

    // Const version of operator[]
    const IndexMap& operator[](size_t index) const {
        assert(index < data.size());
        return data[index];
    }

    // Get matrix size
    size_t size() const { return data.size(); }

    // Iterator types
    using iterator = std::vector<IndexMap>::iterator;
    using const_iterator = std::vector<IndexMap>::const_iterator;

    // Iterator methods
    iterator begin() { return data.begin(); }
    iterator end() { return data.end(); }

    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }

    const_iterator cbegin() const { return data.cbegin(); }
    const_iterator cend() const { return data.cend(); }

    // Frobenius norm (Euclidean norm)
    double norm() const {
        double sum = 0.0;
        for (const auto& row : data) {
            for (const auto& [_, val] : row) {
                sum += val * val;
            }
        }
        return std::sqrt(sum);
    }

    Field diagonal() const {
        Field diagonal(size());
        diagonal.setZero();
        for (size_t i = 0; i < size(); ++i) {
            auto it = data[i].find(i);
            if (it != data[i].end()) {
                diagonal[i] = it->second;
            }
        }
        return diagonal;
    }

    size_t rows() const { return data.size(); }
    size_t cols() const { return data.size(); }
    size_t nnz() const {
        size_t count = 0;
        for (const auto& row : data) {
            count += row.size();
        }
        return count;
    }

    Operator convert_to_eigen() const {
        std::vector<Trip> trips;
        trips.reserve(nnz());
        for (size_t row = 0; row < size(); ++row) {
            // row is IndexMap
            for (auto& [col, value] : data[row]) {
                trips.emplace_back(row, col, value);
            }
        }
        Operator mat(rows(), cols());
        mat.setFromTriplets(trips.begin(), trips.end());
        return mat;
    }

   private:
    std::vector<IndexMap> data;
};
