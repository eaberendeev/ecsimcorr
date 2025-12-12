#ifndef VEC_H_
#define VEC_H_
#include "util.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using  ulog = unsigned int;
// *** vec2 ****
template <typename T>
struct vec2 {
    vec2 (const T u, const T v) : d{u,v} {}
    vec2 (const T a[2]) : d{a[0],a[1]} {}
    vec2 (): d{0,0} {}
    T& operator()(int i) {return d[i];}
    const T& operator()(int i) const {return d[i];}
    vec2<T>& operator=(double s) {d[0]=s;d[1]=s;return (*this);}
    vec2<T>& operator+=(const vec2<T>& o) {d[0]+=o(0);d[1]+=o(1);return(*this);}
    vec2<T>& operator-=(const vec2<T>& o) {d[0]-=o(0);d[1]-=o(1);return(*this);}
    vec2<T> operator/(double s) {vec2<T>o; o(0)=d[0]/s;o(1)=d[1]/s;return o;}
    vec2<T>& operator/=(double s) {d[0]/=s;d[1]/=s;return (*this);}

    //dot product of two vectors
    friend T dot(const vec2<T> &v1, const vec2<T> &v2) {
        T s=0;  for (int i=0;i<2;i++) s+=v1(i)*v2(i);
        return s;   }

    //vector magnitude
    friend T mag(const vec2<T> &v) {return sqrt(dot(v,v));}

    //unit vector
    friend vec2<T> unit(const vec2<T> &v) {return vec2(v)/mag(v);}
    int total_size() const { return d[0] * d[1]; }

    T& x() {return d[0];}
    T x() const {return d[0];}
    T& y() {return d[1];}
    T y() const {return d[1];}
protected:
    T d[2];
};

//vec2-vec2 operations
template<typename T>    //addition of two vec3s
vec2<T> operator+(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)+b(0),a(1)+b(1));   }
template<typename T>    //subtraction of two vec2s
vec2<T> operator-(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)-b(0),a(1)-b(1));   }
template<typename T>    //element-wise multiplication of two vec2s
vec2<T> operator*(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)*b(0),a(1)*b(1));   }
template<typename T>    //element wise division of two vec3s
vec2<T> operator/(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)/b(0),a(1)/b(1));   }

//vec2 - scalar operations
template<typename T>        //scalar multiplication
vec2<T> operator*(const vec2<T> &a, T s) {
    return vec2<T>(a(0)*s, a(1)*s);}
template<typename T>        //scalar multiplication 2
vec2<T> operator*(T s,const vec2<T> &a) {
    return vec2<T>(a(0)*s, a(1)*s);}

template <typename T>
bool operator!=(const vec2<T>& a, const vec2<T>& b) {
    return (a(0) != b(0) || a(1) != b(1));
}
template <typename T>
bool operator==(const vec2<T>& a, const vec2<T>& b) {
    return !(a != b);
}
//output
template<typename T>    //ostream output
std::ostream& operator<<(std::ostream &out, vec2<T>& v) {
    out<<v(0)<<" "<<v(1)<<" 0"; //paraview does not support 2-component arrays
    return out;
}

using double2 = vec2<double>;
using int2 = vec2<int>;
using int2 = vec2<int>;

template <typename T>
struct vec3 {
    vec3 (const T u, const T v, const T w) : d{u,v,w} {}
    vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
    vec3 (): d{0,0,0} {}
    explicit vec3(const T s) : d{s, s, s} {}
    T& operator()(int i) {return d[i];}
    T operator()(int i) const {return d[i];}
    T& operator[](int i) { return d[i]; }
    T operator[](int i) const { return d[i]; }
    vec3<T>& operator=(double s) {d[0]=s;d[1]=s;d[2]=s;return (*this);}
    vec3<T>& operator+=(const vec3<T>& v) {d[0]+=v(0);d[1]+=v(1);d[2]+=v(2);return(*this);}
    vec3<T>& operator-=(const vec3<T>& v) {d[0]-=v(0);d[1]-=v(1);d[2]-=v(2);return(*this);}
    vec3<T> operator/(T s) const {vec3<T> v; v(0) = d[0] / s; v(1) = d[1] / s; v(2) = d[2] / s; return v;}
    vec3<T>& operator/=(double s) {
        d[0] /= s;
        d[1] /= s;
        d[2] /= s;
        return (*this);
    }

    T dot(const vec3& other) const {
        return d[0] * other(0) + d[1] * other(1) + d[2] * other(2);
    }

    double norm() const { return sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]); }

    T squared() const { return dot(*this); }
    vec3<T> cross(const vec3<T>& other) const {
        return vec3<T>{
            +(d[1] * other(2) - d[2] * other(1)),
            -(d[0] * other(2) - d[2] * other(0)),
            +(d[0] * other(1) - d[1] * other(0)),
        };
    }
    vec3<T> parallel_to(const vec3<T>& ref) const {
        return ((*this).dot(ref) * ref) / ref.square();
    }

    vec3<T> transverse_to(const vec3<T>& ref) const {
        return (*this) - parallel_to(ref);
    }

    T length() const { return std::hypot(d[0], d[1], d[2]); }

    vec3<double> normalized() const {
        double l = length();
        if (l > 1.e-16)
            return operator/(l);

        return vec3<double>(0.0, 0.0, 0.0);
    }
    // todo: do it only for T = int
    int elements_product() const { return d[0] * d[1] * d[2]; }
    T& x() {return d[0];}
    T x() const {return d[0];}
    T& y() {return d[1];}
    T y() const {return d[1];}
    T& z() {return d[2];}
    T z() const {return d[2];}
protected:
    T d[3];
};

//vec3-vec3 operations
template<typename T>    //addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)+b(0),a(1)+b(1),a(2)+b(2)); }
template<typename T>    //subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)-b(0),a(1)-b(1),a(2)-b(2)); }
template<typename T>    //element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)*b(0),a(1)*b(1),a(2)*b(2)); }
template<typename T>    //element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)/b(0),a(1)/b(1),a(2)/b(2)); }
template <typename T>   // element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const T& b) {
    return vec3<T>(a(0) / b, a(1) / b, a(2) / b);
}
//vec3 - scalar operations
template<typename T>        //scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
    return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}
template<typename T>        //scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
    return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}

template <typename T>
bool operator!=(const vec3<T>& a, const vec3<T>& b) {
    return (a(0) != b(0) || a(1) != b(1) || a(2) != b(2));
}
template <typename T>
bool operator==(const vec3<T>& a, const vec3<T>& b) {
    return !(a != b);
}
//output
template<typename T>    //ostream output
std::ostream& operator<<(std::ostream &out, const vec3<T>& v) {
    out<<v(0)<<" "<<v(1)<<" "<<v(2);
    return out;
}

//using double3 = vec3<double>;
//using int3 = vec3<int>;

template<typename T>
struct Vector3 {
  static constexpr int dim = 3;
  T data[dim];

  constexpr Vector3() : data{
        static_cast<T>(0),
        static_cast<T>(0),
        static_cast<T>(0),
      } {}

  constexpr Vector3(const T& v)
    : data{v, v, v} {}

  constexpr Vector3(const T& x, const T& y, const T& z)
    : data{x, y, z} {}

  constexpr Vector3(const T v[Vector3::dim])
    : data{v[X], v[Y], v[Z]} {}

  constexpr operator const T*() const { return data; }

  constexpr operator T*() { return data; }

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
  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  Vector3& operator/=(T scalar) {
      if (scalar == static_cast<T>(0)){
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

  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  Vector3<double> operator/(T scalar) const {
      if (scalar == static_cast<T>(0))
          return {};
      return Vector3{
          static_cast<double>(data[X]) / scalar,
          static_cast<double>(data[Y]) / scalar,
          static_cast<double>(data[Z]) / scalar,
      };
  }
  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  Vector3<double> operator/(const Vector3<double>& v) const {
      if (v[X] == static_cast<T>(0) || v[Y] == static_cast<T>(0) ||
          v[Z] == static_cast<T>(0))
          return {};
      return Vector3{
          static_cast<double>(data[X]) / v[X],
          static_cast<double>(data[Y]) / v[Y],
          static_cast<double>(data[Z]) / v[Z],
      };
  }
  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  Vector3<double> normalized() const {
      return operator/(length());
  }

  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  double length() const {
      return std::hypot(data[X], data[Y], data[Z]);
  }
  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  double norm() const {
      return std::hypot(data[X], data[Y], data[Z]);
  }

  template <typename U = T, typename = std::enable_if_t<std::is_integral_v<U>>>
  T elements_product() const {
      return data[X] * data[Y] * data[Z];
  }

  T elements_sum() const { return data[X] + data[Y] + data[Z]; }

  T dot(const Vector3& other) const {
      return                     //
          data[X] * other[X] +   //
          data[Y] * other[Y] +   //
          data[Z] * other[Z];
  }

  T squared() const { return dot(*this); }

  T abs_max() const {
      return std::max(std::abs(data[X]),   //
                      std::max(std::abs(data[Y]), std::abs(data[Z])));
  }
  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  Vector3<double> parallel_to(const Vector3& ref) const {
      return ((*this).dot(ref) * ref) / ref.squared();
  }
  template <typename U = T,
            typename = std::enable_if_t<std::is_floating_point_v<U>>>
  Vector3<double> transverse_to(const Vector3& ref) const {
      return (*this) - parallel_to(ref);
  }

  void swap_order() { std::swap(data[X], data[Z]); }

  Vector3 cross(const Vector3& other) const {
      return Vector3{
          +(data[Y] * other[Z] - data[Z] * other[Y]),
          -(data[X] * other[Z] - data[Z] * other[X]),
          +(data[X] * other[Y] - data[Y] * other[X]),
      };
  }

  friend std::ostream& operator<<(std::ostream& out, const Vector3& vector) {
      out << std::to_string(vector[X]) << " ";
      out << std::to_string(vector[Y]) << " ";
      out << std::to_string(vector[Z]) << " ";
      return out;
  }
};

//using Vector3R = Vector3<double>;
using double3 = Vector3<double>;
//using Vector3I = Vector3<int>;
using int3 = Vector3<int>;

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

#endif // VEC_H
