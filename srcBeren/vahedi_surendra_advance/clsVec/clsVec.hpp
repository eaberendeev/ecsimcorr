// Реализация класса вектора
#ifndef CLS_VEC_H
#define CLS_VEC_H

#include <cmath>
#include <iostream>

using namespace std;

// Класс вектора
template <typename T>
class clsVec {
    
private:

    T _coord[3];

public:

    clsVec();
    clsVec(T, T, T);
    clsVec(T, size_t);

    T x() const;
    T y() const;
    T z() const;

    T& x();
    T& y();
    T& z();

    void x(T);
    void y(T);
    void z(T);

    clsVec<T> operator+(const clsVec<T>&) const;
    clsVec<T> operator-(const clsVec<T>&) const;
    clsVec<T> operator+() const;
    clsVec<T> operator-() const;
    T& operator[](const int&);
};

extern template class clsVec<size_t>;
extern template class clsVec<float>;
extern template class clsVec<double>;

template <typename T>
T sp(const clsVec<T>& fir_arg, const clsVec<T>& sec_arg) {
    return fir_arg.x() * sec_arg.x() +
        fir_arg.y() * sec_arg.y() +
        fir_arg.z() * sec_arg.z();
            
}

template <typename T>
clsVec<T> vp(const clsVec<T>& fir_arg, const clsVec<T>& sec_arg) {
    return {fir_arg.y() * sec_arg.z() - fir_arg.z() * sec_arg.y(),
            fir_arg.z() * sec_arg.x() - fir_arg.x() * sec_arg.z(),
            fir_arg.x() * sec_arg.y() - fir_arg.y() * sec_arg.x()};
}

template <typename T>
T norm(const clsVec<T>& fir_arg) {
    return sp(fir_arg, fir_arg);
}

template <typename T>
clsVec<T> RotationZ(const clsVec<T>& arg, T phi) {
    return {arg.x()*cos(phi)-arg.y()*sin(phi), arg.x()*sin(phi)+arg.y()*cos(phi), arg.z()};
}

template <typename T>
clsVec<T> RotationY(const clsVec<T>& arg, T phi) {
    return {arg.x()*cos(phi)+arg.z()*sin(phi), arg.y(), -arg.x()*sin(phi)+arg.z()*cos(phi)};
}

template <typename T1, typename T2>   
clsVec<T2> operator*(T1 a, const clsVec<T2>& arg) {
    return {a * arg.x(), a * arg.y(), a * arg.z()};
} 

template <typename T1, typename T2> 
clsVec<T2> operator*(const clsVec<T2>& arg, T1 a) {
    return {a * arg.x(), a * arg.y(), a * arg.z()};
}

template <typename T1, typename T2> 
clsVec<T2> operator/(const clsVec<T2>& arg, T1 a) {
    return {arg.x()/a, arg.y() / a, arg.z() / a};
}

#endif