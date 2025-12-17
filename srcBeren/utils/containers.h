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

#include "util.h"
#include "vector2.h"
#include "vector3.h"
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
    Array3D(const Vector3I& nn) { resize(nn[0], nn[1], nn[2]); }
    Array3D(Array3D&& other) : _data{other._data}, _size{other._size} {}

    Array3D() { _size = Vector3I(0, 0, 0); }

    void resize(int n1, int n2, int n3) {
        _data.resize(n1 * n2 * n3);
        _size = Vector3I(n1, n2, n3);
    }
    void resize(const Vector3I& nn) {
        _data.resize(nn.elements_product());
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
        assert(i < _size[0] && j < _size[1] && k < _size[2]);
        assert(i >= 0 && j >= 0 && k >= 0);
        return true;
    }

    T& operator()(int i, int j, int k) {
        assert(checkIndex(i, j, k));
        return _data[i * _size[1] * _size[2] + j * _size[2] + k];
    }

    const T& operator()(int i, int j, int k) const {
        assert(checkIndex(i, j, k));
        return _data[i * _size[1] * _size[2] + j * _size[2] + k];
    }
    T& operator()(int i) { return _data[i]; }
    const T& operator()(int i) const { return _data[i]; }
    std::vector<T>& data() { return _data; }
    const std::vector<T>& data() const { return _data; }

    Vector3I size() const { return _size; }
    int capacity() const { return _size.elements_product(); }

   private:
    std::vector<T> _data;
    Vector3I _size;
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
        assert(i < _size[0] && j < _size[1]);
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
    Field3d(const Vector3I& nn, int d) { resize(nn, d); }
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
        other._size = Vector3I(0, 0, 0);
        other._nd = 0;
    }
    Field3d() { _size = Vector3I(0, 0, 0); }
    // Add static Zero constructor
    static Field3d Zero(const Vector3I& size, int nd) {
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
        _size = Vector3I(n1, n2, n3);
        _nd = d;
    }
    void resize(const Vector3I& nn, int d) {
        _data.resize(nn.elements_product() * d);
        _size = nn;
        _nd = d;
    }
    void resize(int n){
        _data.resize(n);
        _size = Vector3I(n, 1, 1);
        _nd = 1;
    }
    void setZero() {
#pragma omp parallel for simd
        for (auto i = 0; i < capacity(); ++i) {
            _data[i] = 0.;
        }
    }
    bool checkIndex(int i, int j, int k, int d) const {
        bool true1 = i < _size[0] && j < _size[1] && k < _size[2] && d < _nd;
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
            other._size = Vector3I(0, 0, 0);
            other._nd = 0;
        }
        return *this;
    }
    double& operator()(int i, int j, int k, int d) {
        assert(checkIndex(i, j, k, d));
        return _data[d + _nd * (i * _size[1] * _size[2] + j * _size[2] + k)];
    }
    const double& operator()(int i, int j, int k, int d) const {
        assert(checkIndex(i, j, k, d));
        return _data[d + _nd * (i * _size[1] * _size[2] + j * _size[2] + k)];
    }
    double& operator()(int i) { return _data[i]; }
    const double& operator()(int i) const { return _data[i]; }
    double& operator[](int i) { return _data[i]; }
    const double& operator[](int i) const { return _data[i]; }
    int nd() const { return _nd; }
    Vector3I sizes() const { return _size; }
    int capacity() const { return _nd * _size.elements_product(); }
    size_t size() const { return (size_t) _nd * _size.elements_product(); }
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
    double squared() const {
        double result = 0.0;
#pragma omp parallel for simd reduction(+ : result)
     for (int i = 0; i < capacity(); ++i) {
        result += _data[i] * _data[i];
    }
    return result;
    }

    double norm() const { return std::sqrt(squared()); }

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
    Vector3I _size;
    int _nd;
};

#endif
