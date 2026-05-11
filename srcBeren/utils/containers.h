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

#include "vector2.h"
#include "vector3.h"
#include "indexing.h"
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

[[noreturn]] inline void fatalDimensionMismatch(const char* context) {
    std::cerr << "Fatal error: array dimension mismatch in " << context
              << std::endl;
    std::abort();
}
template <typename T>
class Array3D {
   public:
    Array3D() : size_() {}
    Array3D(int n1, int n2, int n3) { resize(n1, n2, n3); }
    explicit Array3D(const Vector3I& nn) { resize(nn); }

    Array3D(const Array3D&) = default;

    Array3D(Array3D&& other) noexcept
        : data_(std::move(other.data_)), size_(other.size_) {
        other.size_ = {};
    }

    Array3D& operator=(const Array3D& other) {
        if (this == &other)
            return *this;
        if (size_ != other.size_) {
            fatalDimensionMismatch("Array3D copy assignment");
        }
        data_ = other.data_;
        return *this;
    }

    Array3D& operator=(Array3D&& other) {
        if (this == &other)
            return *this;
        if (size_ != other.size_) {
            fatalDimensionMismatch("Array3D move assignment");
        }
        data_ = std::move(other.data_);
        size_ = other.size_;
        other.size_ = {};
        return *this;
    }

    void resize(int n1, int n2, int n3) {
        data_.resize(static_cast<std::size_t>(n1) * n2 * n3);
        size_ = Vector3I(n1, n2, n3);
    }

    void resize(const Vector3I& nn) {
        data_.resize(nn.elements_product());
        size_ = nn;
    }

    void setZero() { std::fill(data_.begin(), data_.end(), T(0)); }

    T& operator()(int i, int j, int k) {
        assert(checkIndex(i, j, k));
        return data_[static_cast<std::size_t>(i) * size_[1] * size_[2] +
                     static_cast<std::size_t>(j) * size_[2] + k];
    }

    const T& operator()(int i, int j, int k) const {
        assert(checkIndex(i, j, k));
        return data_[static_cast<std::size_t>(i) * size_[1] * size_[2] +
                     static_cast<std::size_t>(j) * size_[2] + k];
    }

    T& operator()(int i) { return data_[i]; }
    const T& operator()(int i) const { return data_[i]; }

    std::vector<T>& data() { return data_; }
    const std::vector<T>& data() const { return data_; }

    Vector3I size() const { return size_; }
    int capacity() const { return size_.elements_product(); }

   private:
    std::vector<T> data_;
    Vector3I size_;

    bool checkIndex(int i, int j, int k) const {
        assert(i >= 0 && i < size_[0]);
        assert(j >= 0 && j < size_[1]);
        assert(k >= 0 && k < size_[2]);
        return true;
    }
};

template <typename T>
class Array2D {
   public:
    Array2D() : size_() {}
    Array2D(int n1, int n2) { resize(n1, n2); }
    explicit Array2D(const int2& nn) { resize(nn); }

    Array2D(const Array2D&) = default;

    Array2D(Array2D&& other) noexcept
        : data_(std::move(other.data_)), size_(other.size_) {
        other.size_ = {};
    }

    Array2D& operator=(const Array2D& other) {
        if (this == &other)
            return *this;
        if (size_ != other.size_) {
            fatalDimensionMismatch("Array2D copy assignment");
        }
        data_ = other.data_;
        return *this;
    }

    Array2D& operator=(Array2D&& other) {
        if (this == &other)
            return *this;
        if (size_ != other.size_) {
            fatalDimensionMismatch("Array2D move assignment");
        }
        data_ = std::move(other.data_);
        size_ = other.size_;
        other.size_ = {};
        return *this;
    }

    void resize(int n1, int n2) {
        data_.resize(static_cast<std::size_t>(n1) * n2);
        size_ = int2(n1, n2);
    }

    void resize(const int2& nn) {
        data_.resize(nn.total_size());
        size_ = nn;
    }

    void setZero() { std::fill(data_.begin(), data_.end(), T(0)); }

    T& operator()(int i, int j) {
        assert(checkIndex(i, j));
        return data_[static_cast<std::size_t>(i) * size_[1] + j];
    }

    const T& operator()(int i, int j) const {
        assert(checkIndex(i, j));
        return data_[static_cast<std::size_t>(i) * size_[1] + j];
    }

    T& operator()(int i) { return data_[i]; }
    const T& operator()(int i) const { return data_[i]; }

    std::vector<T>& data() { return data_; }
    const std::vector<T>& data() const { return data_; }

    int2 size() const { return size_; }
    int capacity() const { return size_.total_size(); }

   private:
    std::vector<T> data_;
    int2 size_;

    bool checkIndex(int i, int j) const {
        assert(i >= 0 && i < size_[0]);
        assert(j >= 0 && j < size_[1]);
        return true;
    }
};
class Field3d {
   public:
    Field3d(const int n) { resize(n); }
    Field3d(int n1, int n2, int n3, int d) { resize(n1, n2, n3, d); }
    Field3d(const Vector3I& nn, int d) { resize(nn, d); }
    Field3d(const Field3d& other)
        : data_(other.data_), size_(other.size_), nd_(other.nd_) {}

    Field3d(Field3d&& other) noexcept
        : data_(std::move(other.data_)), size_(other.size_), nd_(other.nd_) {
        other.size_ = Vector3I(0, 0, 0);
        other.nd_ = 0;
    }

    Field3d() : size_(0, 0, 0), nd_(0) {}

    static Field3d Zero(const Vector3I& size, int nd) {
        Field3d result(size, nd);
        result.setZero();
        return result;
    }
    static Field3d Zero(int n) {
        Field3d result(n);
        result.setZero();
        return result;
    }
    void resize(int n1, int n2, int n3, int d) {
        data_.resize(n1 * n2 * n3 * d);
        size_ = Vector3I(n1, n2, n3);
        nd_ = d;
    }
    void resize(const Vector3I& nn, int d) {
        data_.resize(nn.elements_product() * d);
        size_ = nn;
        nd_ = d;
    }
    void resize(int n) {
        data_.resize(n);
        size_ = Vector3I(n, 1, 1);
        nd_ = 1;
    }
    void setZero() {
        data_.setZero();
    }

    bool checkIndex(int i, int j, int k, int d) const {
        return i >= 0 && i < size_[0] && j >= 0 && j < size_[1] && k >= 0 &&
               k < size_[2] && d >= 0 && d < nd_;
    }
    Field3d& operator=(const Field3d& other) {
        if (this == &other)
            return *this;
        if (size_ != other.size_ || nd_ != other.nd_)
            fatalDimensionMismatch("copy assignment");
        data_ = other.data_;
        return *this;
    }

    Field3d& operator=(Field3d&& other) noexcept {
        if (this != &other) {
            data_ = std::move(other.data_);
            size_ = other.size_;
            nd_ = other.nd_;
            other.size_ = Vector3I(0, 0, 0);
            other.nd_ = 0;
        }
        return *this;
    }

    double& operator()(int i, int j, int k, int d) {
        assert(checkIndex(i, j, k, d));
        return data_[d + nd_ * (i * size_[1] * size_[2] + j * size_[2] + k)];
    }
    const double& operator()(int i, int j, int k, int d) const {
        assert(checkIndex(i, j, k, d));
        return data_[d + nd_ * (i * size_[1] * size_[2] + j * size_[2] + k)];
    }

    double& operator()(int i) { return data_[i]; }
    const double& operator()(int i) const { return data_[i]; }
    double& operator[](int i) { return data_[i]; }
    const double& operator[](int i) const { return data_[i]; }

    int nd() const { return nd_; }
    Vector3I sizes() const { return size_; }
    int capacity() const { return nd_ * size_.elements_product(); }
    size_t size() const {
        return static_cast<size_t>(nd_) * size_.elements_product();
    }

    Eigen::VectorXd& data() { return data_; }
    const Eigen::VectorXd& data() const { return data_; }

    Field3d& operator*=(const double alpha) {
        data_ *= alpha;
        return *this;
    }

    Field3d& operator+=(const Field3d& other) {
        if (capacity() != other.capacity())
            fatalDimensionMismatch("operator+=");
        data_ += other.data_;
        return *this;
    }

    Field3d& operator-=(const Field3d& other) {
        if (capacity() != other.capacity())
            fatalDimensionMismatch("operator-=");
        data_ -= other.data_;
        return *this;
    }

    double dot(const Field3d& other) const {
        if (capacity() != other.capacity())
            fatalDimensionMismatch("dot");
        return data_.dot(other.data_);
    }

    double squared() const { return data_.squaredNorm(); }

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
        Field3d res(field.sizes(), field.nd());
#pragma omp parallel for schedule(dynamic, 256)
        for (int row = 0; row < rows; ++row) {
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

    using iterator = Eigen::VectorXd::iterator;
    using const_iterator = Eigen::VectorXd::const_iterator;

    iterator begin() { return data_.begin(); }
    iterator end() { return data_.end(); }
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator cend() const { return data_.cend(); }

   private:
    Eigen::VectorXd data_;
    Vector3I size_;
    int nd_;
};

static inline double dot_product_sum(const Field3d& f, const Field3d& g,
                                     const IndexRange& range) {
    double accumulator = 0;

    for (auto i = range.start.x(); i < range.end.x(); ++i) {
        for (auto j = range.start.y(); j < range.end.y(); ++j) {
            for (auto k = range.start.z(); k < range.end.z(); ++k) {
                Vector3R v1 =
                    Vector3R(f(i, j, k, 0), f(i, j, k, 1), f(i, j, k, 2));
                Vector3R v2 =
                    Vector3R(g(i, j, k, 0), g(i, j, k, 1), g(i, j, k, 2));
                accumulator += v1.dot(v2);
            }
        }
    }
    return accumulator;
}

#endif
