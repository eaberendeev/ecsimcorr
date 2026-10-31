#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>

#include "core/types.h"

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

    constexpr Vector3(const T v[Vector3::dim]) : data{v[0], v[1], v[2]} {
    }

    constexpr operator const T*() const {
        return data;
    }

    constexpr operator T*() {
        return data;
    }

    // clang-format off: access specifiers
  constexpr T& x() { return data[0]; }
  constexpr T& y() { return data[1]; }
  constexpr T& z() { return data[2]; }
  constexpr const T& x() const { return data[0]; }
  constexpr const T& y() const { return data[1]; }
  constexpr const T& z() const { return data[2]; }

  constexpr T& operator[](int i) { return data[i]; }
  constexpr const T& operator[](int i) const { return data[i]; }
    // clang-format on

    bool operator==(const Vector3& other) const {
        return                       //
            data[0] == other[0] &&   //
            data[1] == other[1] &&   //
            data[2] == other[2];     //
    }
    bool operator!=(const Vector3& other) const {
        return !(*this == other);
    }

    Vector3& operator+=(const Vector3& other) {
        data[0] += other[0];
        data[1] += other[1];
        data[2] += other[2];
        return *this;
    }

    Vector3& operator-=(const Vector3& other) {
        data[0] -= other[0];
        data[1] -= other[1];
        data[2] -= other[2];
        return *this;
    }

    Vector3& operator*=(T scalar) {
        data[0] *= scalar;
        data[1] *= scalar;
        data[2] *= scalar;
        return *this;
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3& operator/=(T scalar) {
        if (scalar == static_cast<T>(0)) {
            data[0] = data[1] = data[2] = 0;
            return *this;
        }
        data[0] /= scalar;
        data[1] /= scalar;
        data[2] /= scalar;
        return *this;
    }

    Vector3 operator+(const Vector3& other) const {
        return Vector3{
            data[0] + other[0],
            data[1] + other[1],
            data[2] + other[2],
        };
    }

    Vector3 operator-(const Vector3& other) const {
        return Vector3{
            data[0] - other[0],
            data[1] - other[1],
            data[2] - other[2],
        };
    }

    Vector3 elementwise_product(const Vector3& other) const {
        return Vector3{
            data[0] * other[0],
            data[1] * other[1],
            data[2] * other[2],
        };
    }

    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> operator/(T scalar) const {
        if (scalar == static_cast<T>(0))
            return {};
        return Vector3{
            static_cast<double>(data[0]) / scalar,
            static_cast<double>(data[1]) / scalar,
            static_cast<double>(data[2]) / scalar,
        };
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> operator/(const Vector3<double>& v) const {
        if (v[0] == static_cast<T>(0) || v[1] == static_cast<T>(0) || v[2] == static_cast<T>(0))
            return {};
        return Vector3{
            static_cast<double>(data[0]) / v[0],
            static_cast<double>(data[1]) / v[1],
            static_cast<double>(data[2]) / v[2],
        };
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> normalized() const {
        return operator/(length());
    }

    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    double length() const {
        return std::hypot(data[0], data[1], data[2]);
    }
    //   template <typename U = T,
    //             typename = std::enable_if_t<std::is_floating_point_v<U>>>
    double norm() const {
        return std::hypot(data[0], data[1], data[2]);
    }

    T elements_product() const {
        return data[0] * data[1] * data[2];
    }

    T elements_sum() const {
        return data[0] + data[1] + data[2];
    }

    T dot(const Vector3& other) const {
        return                     //
            data[0] * other[0] +   //
            data[1] * other[1] +   //
            data[2] * other[2];
    }

    T squaredNorm() const {
        return dot(*this);
    }

    T abs_max() const {
        return std::max(std::abs(data[0]),   //
                        std::max(std::abs(data[1]), std::abs(data[2])));
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> parallel_to(const Vector3& ref) const {
        return ((*this).dot(ref) * ref) / ref.squaredNorm();
    }
    template <typename U = T, typename = std::enable_if_t<std::is_floating_point_v<U>>>
    Vector3<double> transverse_to(const Vector3& ref) const {
        return (*this) - parallel_to(ref);
    }

    void swap_order() {
        std::swap(data[0], data[2]);
    }

    Vector3 cross(const Vector3& other) const {
        return Vector3{
            +(data[1] * other[2] - data[2] * other[1]),
            -(data[0] * other[2] - data[2] * other[0]),
            +(data[0] * other[1] - data[1] * other[0]),
        };
    }
    auto split() const {
        return std::array<T, 3>{data[0], data[1], data[2]};
    }

    friend std::ostream& operator<<(std::ostream& out, const Vector3& vector) {
        out << std::to_string(vector[0]) << " ";
        out << std::to_string(vector[1]) << " ";
        out << std::to_string(vector[2]) << " ";
        return out;
    }
};

using Vector3R = Vector3<double>;
using Vector3I = Vector3<int>;

template <typename T>
Vector3<T> operator*(const Vector3<T>& vector, T scalar) {
    return Vector3{
        vector[0] * scalar,
        vector[1] * scalar,
        vector[2] * scalar,
    };
}

template <typename T>
Vector3<T> operator*(T scalar, const Vector3<T>& vector) {
    return vector * scalar;
}

template <typename T>
Vector3<T> min(const Vector3<T>& lhs, const Vector3<T>& rhs) {
    return Vector3{
        std::min(lhs[0], rhs[0]),
        std::min(lhs[1], rhs[1]),
        std::min(lhs[2], rhs[2]),
    };
}

template <typename T>
Vector3<T> max(const Vector3<T>& lhs, const Vector3<T>& rhs) {
    return Vector3{
        std::max(lhs[0], rhs[0]),
        std::max(lhs[1], rhs[1]),
        std::max(lhs[2], rhs[2]),
    };
}
