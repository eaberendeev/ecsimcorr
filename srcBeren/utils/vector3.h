#pragma once

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>
#include <string>
#include <vector>

#include "types.h"

template <typename T>
struct Vector3 {
    static constexpr int dim = 3;
    T data[dim];

    constexpr Vector3()
        : data{
              static_cast<T>(0),
              static_cast<T>(0),
              static_cast<T>(0),
          } {
    }

    constexpr Vector3(const T& v) : data{v, v, v} {
    }

    constexpr Vector3(const T& x, const T& y, const T& z) : data{x, y, z} {
    }

    constexpr Vector3(const T v[Vector3::dim]) : data{v[X], v[Y], v[Z]} {
    }

    constexpr operator const T*() const {
        return data;
    }

    constexpr operator T*() {
        return data;
    }

    // clang-format off: access specifiers
  constexpr T& x() { return data[X]; }
  constexpr T& y() { return data[Y]; }
  constexpr T& z() { return data[Z]; }
  constexpr const T& x() const { return data[X]; }
  constexpr const T& y() const { return data[Y]; }
  constexpr const T& z() const { return data[Z]; }

  constexpr T& operator[](int i) { return data[i]; }
  constexpr const T& operator[](int i) const { return data[i]; }
    // clang-format on

    bool operator==(const Vector3& other) const {
        return                       //
            data[X] == other[X] &&   //
            data[Y] == other[Y] &&   //
            data[Z] == other[Z];     //
    }
    bool operator!=(const Vector3& other) const {
        return !(*this == other);
    }

    Vector3& operator+=(const Vector3& other) {
        data[X] += other[X];
        data[Y] += other[Y];
        data[Z] += other[Z];
        return *this;
    }

    Vector3& operator-=(const Vector3& other) {
        data[X] -= other[X];
        data[Y] -= other[Y];
        data[Z] -= other[Z];
        return *this;
    }

    Vector3& operator*=(T scalar) {
        data[X] *= scalar;
        data[Y] *= scalar;
        data[Z] *= scalar;
        return *this;
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3& operator/=(T scalar) {
        if (scalar == static_cast<T>(0)) {
            data[X] = data[Y] = data[Z] = 0;
            return *this;
        }
        data[X] /= scalar;
        data[Y] /= scalar;
        data[Z] /= scalar;
        return *this;
    }

    Vector3 operator+(const Vector3& other) const {
        return Vector3{
            data[X] + other[X],
            data[Y] + other[Y],
            data[Z] + other[Z],
        };
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3{
            data[X] - other[X],
            data[Y] - other[Y],
            data[Z] - other[Z],
        };
    }

    Vector3 elementwise_product(const Vector3& other) const {
        return Vector3{
            data[X] * other[X],
            data[Y] * other[Y],
            data[Z] * other[Z],
        };
    }

    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> operator/(T scalar) const {
        if (scalar == static_cast<T>(0))
            return {};
        return Vector3{
            static_cast<double>(data[X]) / scalar,
            static_cast<double>(data[Y]) / scalar,
            static_cast<double>(data[Z]) / scalar,
        };
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> operator/(const Vector3<double>& v) const {
        if (v[X] == static_cast<T>(0) || v[Y] == static_cast<T>(0) || v[Z] == static_cast<T>(0))
            return {};
        return Vector3{
            static_cast<double>(data[X]) / v[X],
            static_cast<double>(data[Y]) / v[Y],
            static_cast<double>(data[Z]) / v[Z],
        };
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> normalized() const {
        return operator/(length());
    }

    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    double length() const {
        return std::hypot(data[X], data[Y], data[Z]);
    }
    //   template <typename U = T,
    //             typename = std::enable_if_t<std::is_floating_point_v<U>>>
    double norm() const {
        return std::hypot(data[X], data[Y], data[Z]);
    }

    T elements_product() const {
        return data[X] * data[Y] * data[Z];
    }

    T elements_sum() const {
        return data[X] + data[Y] + data[Z];
    }

    T dot(const Vector3& other) const {
        return                     //
            data[X] * other[X] +   //
            data[Y] * other[Y] +   //
            data[Z] * other[Z];
    }

    T squared() const {
        return dot(*this);
    }

    T abs_max() const {
        return std::max(std::abs(data[X]),   //
                        std::max(std::abs(data[Y]), std::abs(data[Z])));
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> parallel_to(const Vector3& ref) const {
        return ((*this).dot(ref) * ref) / ref.squared();
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> transverse_to(const Vector3& ref) const {
        return (*this) - parallel_to(ref);
    }

    void swap_order() {
        std::swap(data[X], data[Z]);
    }

    Vector3 cross(const Vector3& other) const {
        return Vector3{
            +(data[Y] * other[Z] - data[Z] * other[Y]),
            -(data[X] * other[Z] - data[Z] * other[X]),
            +(data[X] * other[Y] - data[Y] * other[X]),
        };
    }
    auto split() const {
        return std::array<T, 3>{data[X], data[Y], data[Z]};
    }

    friend std::ostream& operator<<(std::ostream& out, const Vector3& vector) {
        out << std::to_string(vector[X]) << " ";
        out << std::to_string(vector[Y]) << " ";
        out << std::to_string(vector[Z]) << " ";
        return out;
    }
};

using Vector3R = Vector3<double>;
using Vector3I = Vector3<int>;

template <typename T>
Vector3<T> operator*(const Vector3<T>& vector, T scalar) {
    return Vector3{
        vector[X] * scalar,
        vector[Y] * scalar,
        vector[Z] * scalar,
    };
}

template <typename T>
Vector3<T> operator*(T scalar, const Vector3<T>& vector) {
    return vector * scalar;
}

template <typename T>
Vector3<T> min(const Vector3<T>& lhs, const Vector3<T>& rhs) {
    return Vector3{
        std::min(lhs[X], rhs[X]),
        std::min(lhs[Y], rhs[Y]),
        std::min(lhs[Z], rhs[Z]),
    };
}

template <typename T>
Vector3<T> max(const Vector3<T>& lhs, const Vector3<T>& rhs) {
    return Vector3{
        std::max(lhs[X], rhs[X]),
        std::max(lhs[Y], rhs[Y]),
        std::max(lhs[Z], rhs[Z]),
    };
}
