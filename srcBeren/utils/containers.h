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
    void set_zero() {
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
    void set_zero() {
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
    Field3d(int n1, int n2, int n3, int d) { resize(n1, n2, n3, d); }
    Field3d(const int3& nn, int d) { resize(nn, d); }
    Field3d(const Field3d& other) : _size{other._size}, _nd{other._nd} {
        resize(_size, _nd);
        _data = other._data;
    }
    Field3d() { _size = int3(0, 0, 0); }

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
    void set_zero() {
        for (auto i = 0; i < capacity(); ++i) {
            _data[i] = 0.;
        }
    }
    bool checkIndex(int i, int j, int k, int d) const {
        assert(i < _size(0) && j < _size(1) && k < _size(2) && d < _nd);
        assert(i >= 0 && j >= 0 && k >= 0 && d >= 0);
        return true;
    }
    Field3d& operator=(Field3d& other) {
        // Self-assignment check
        if (this == &other)
            return *this;
        // Checking the conformity of dimensions
        assert(_size == other._size && _nd == other._nd);
        _data = other._data;
        return (*this);
    }
    Field3d& operator=(Field3d&& other) {
        // Self-assignment check
        if (this == &other)
            return *this;
        // Checking the conformity of dimensions
        assert(_size == other._size && _nd == other._nd);
        _data = other._data;
        _size = other._size;
        _nd = other._nd;
        return (*this);
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
    int nd() const { return _nd; }
    int3 size() const { return _size; }
    int capacity() const { return _nd * _size.total_size(); }
    Eigen::VectorXd& data() { return _data; }
    const Eigen::VectorXd& data() const { return _data; }

   private:
    Eigen::VectorXd _data;
    int3 _size;
    int _nd;
};

#endif
