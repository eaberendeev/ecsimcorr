#ifndef VEC_H_
#define VEC_H_
#include "defines.h"
#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

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
    vec2<T>& operator+=(vec2<T> o) {d[0]+=o(0);d[1]+=o(1);return(*this);}
    vec2<T>& operator-=(vec2<T> o) {d[0]-=o(0);d[1]-=o(1);return(*this);}
    vec2<T> operator/(double s) {vec2<T>o; o(0)=d[0]/s;o(1)=d[1]/s;return o;}
    vec2<T> operator/=(double s) {d[0]/=s;d[1]/=s;return (*this);}

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

    T& operator()(int i) {return d[i];}
    T operator()(int i) const {return d[i];}
    vec3<T>& operator=(double s) {d[0]=s;d[1]=s;d[2]=s;return (*this);}
    vec3<T>& operator+=(vec3<T> v) {d[0]+=v(0);d[1]+=v(1);d[2]+=v(2);return(*this);}
    vec3<T>& operator-=(vec3<T> v) {d[0]-=v(0);d[1]-=v(1);d[2]-=v(2);return(*this);}
    vec3<T> operator/(double s) {vec3<T> v; v(0) = d[0] / s; v(1) = d[1] / s; v(2) = d[2] / s; return v;}
    vec3<T> operator/=(double s) {d[0]/=s;d[1]/=s;d[2]/=s;return (*this);}

    //dot product of two vectors
    friend T dot(const vec3<T> &v1, const vec3<T> &v2) {
        T s=0;  for (int i=0;i<3;i++) s+=v1(i)*v2(i);
        return s;   }

    //vector magnitude
    friend T mag(const vec3<T> &v) {return sqrt(dot(v,v));}

    //unit vector
    friend vec3<T> unit(const vec3<T> &v) {if (mag(v)>0) return vec3(v)/mag(v); else return {0,0,0};}

    //cross product
    friend vec3<T> cross(const vec3<T> &a, const vec3<T> &b) {
        return {a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)};
    }
    int total_size() const { return d[0] * d[1] * d[2]; }
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

using double3 = vec3<double>;
using int3 = vec3<int>;

#endif // VEC_H
