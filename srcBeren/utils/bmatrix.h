#pragma once

#include <assert.h>

#include <unordered_map>
#include <vector>

typedef std::unordered_map<int, double> IndexMap;
typedef std::vector<IndexMap> BMatrix;

// Matrix multiplication:
BMatrix operator*(const BMatrix& a, const BMatrix& b);

// Matrix addition:
BMatrix operator+(const BMatrix& a, const BMatrix& b);

// Matrix subtraction:
BMatrix operator-(const BMatrix& a, const BMatrix& b);

// Matrix-vector multiplication:
std::vector<double> operator*(const BMatrix& mat,
                              const std::vector<double>& vec);
