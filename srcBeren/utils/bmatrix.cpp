#include <assert.h>

#include "bmatrix.h"
#include <vector>

// Matrix multiplication:
BMatrix operator*(const BMatrix& a, const BMatrix& b) {
    BMatrix result(a.size());

    for (size_t i = 0; i < a.size(); i++) {
        for (const auto& [k1, v1] : a[i]) {
            for (const auto& [k2, v2] : b[k1]) {
                result[i][k2] += v1 * v2;
            }
        }
    }
    return result;
}

// Matrix addition:
BMatrix operator+(const BMatrix& a, const BMatrix& b) {
    assert(a.size() == b.size());

    BMatrix result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i].reserve(std::max(a[i].size(), b[i].size()));
    }

    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i];
        for (const auto& [k, v] : b[i]) {
            result[i][k] += v;
        }
    }
    return result;
}

// Matrix subtraction:
BMatrix operator-(const BMatrix& a, const BMatrix& b) {
    BMatrix result(a.size());

    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i];
        for (const auto& [k, v] : b[i]) {
            result[i][k] -= v;
        }
    }
    return result;
}

// Matrix-vector multiplication:
std::vector<double> operator*(const BMatrix& mat,
                              const std::vector<double>& vec) {
    std::vector<double> result(mat.size(), 0.0);

    for (size_t i = 0; i < mat.size(); i++) {
        for (const auto& [k, v] : mat[i]) {
            result[i] += v * vec[k];
        }
    }
    return result;
}
