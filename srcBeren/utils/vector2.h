#pragma once

#include <iostream>
#include <sstream>

// *** Vector2 ****
template <typename T>
struct Vector2 {
    Vector2(const T u, const T v) : d{u, v} {}
    Vector2(const T a[2]) : d{a[0], a[1]} {}
    Vector2() : d{0, 0} {}
    T& operator()(int i) { return d[i]; }
    const T& operator()(int i) const { return d[i]; }
    Vector2<T>& operator=(double s) {
        d[0] = s;
        d[1] = s;
        return (*this);
    }
    Vector2<T>& operator+=(const Vector2<T>& o) {
        d[0] += o(0);
        d[1] += o(1);
        return (*this);
    }
    Vector2<T>& operator-=(const Vector2<T>& o) {
        d[0] -= o(0);
        d[1] -= o(1);
        return (*this);
    }
    Vector2<T> operator/(double s) {
        Vector2<T> o;
        o(0) = d[0] / s;
        o(1) = d[1] / s;
        return o;
    }
    Vector2<T>& operator/=(double s) {
        d[0] /= s;
        d[1] /= s;
        return (*this);
    }

    // dot product of two vectors
    friend T dot(const Vector2<T>& v1, const Vector2<T>& v2) {
        T s = 0;
        for (int i = 0; i < 2; i++) s += v1(i) * v2(i);
        return s;
    }

    // vector magnitude
    friend T mag(const Vector2<T>& v) { return sqrt(dot(v, v)); }

    // unit vector
    friend Vector2<T> unit(const Vector2<T>& v) { return Vector2(v) / mag(v); }
    int total_size() const { return d[0] * d[1]; }

    T& x() { return d[0]; }
    T x() const { return d[0]; }
    T& y() { return d[1]; }
    T y() const { return d[1]; }

   protected:
    T d[2];
};

// Vector2-Vector2 operations
template <typename T>   // addition of two vec3s
Vector2<T> operator+(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>(a(0) + b(0), a(1) + b(1));
}
template <typename T>   // subtraction of two vec2s
Vector2<T> operator-(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>(a(0) - b(0), a(1) - b(1));
}
template <typename T>   // element-wise multiplication of two vec2s
Vector2<T> operator*(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>(a(0) * b(0), a(1) * b(1));
}
template <typename T>   // element wise division of two vec3s
Vector2<T> operator/(const Vector2<T>& a, const Vector2<T>& b) {
    return Vector2<T>(a(0) / b(0), a(1) / b(1));
}

// Vector2 - scalar operations
template <typename T>   // scalar multiplication
Vector2<T> operator*(const Vector2<T>& a, T s) {
    return Vector2<T>(a(0) * s, a(1) * s);
}
template <typename T>   // scalar multiplication 2
Vector2<T> operator*(T s, const Vector2<T>& a) {
    return Vector2<T>(a(0) * s, a(1) * s);
}

template <typename T>
bool operator!=(const Vector2<T>& a, const Vector2<T>& b) {
    return (a(0) != b(0) || a(1) != b(1));
}
template <typename T>
bool operator==(const Vector2<T>& a, const Vector2<T>& b) {
    return !(a != b);
}
// output
template <typename T>   // ostream output
std::ostream& operator<<(std::ostream& out, Vector2<T>& v) {
    out << v(0) << " " << v(1)
        << " 0";   // paraview does not support 2-component arrays
    return out;
}

using double2 = Vector2<double>;
using int2 = Vector2<int>;
