#include <clsVec.hpp>

#include <cmath>
#include <iostream>

using namespace std;

template <typename T>
clsVec<T>::clsVec() {
    _coord[0] = T(0);
    _coord[1] = T(0);
    _coord[2] = T(0);
}

template <typename T>
clsVec<T>::clsVec(T xi, T yi, T zi) {
    _coord[0] = xi;
    _coord[1] = yi;
    _coord[2] = zi;
}

// template <typename T>
// clsVec<T>::clsVec(T coord, size_t axis) {
//     coord[0] = T(0);
//     coord[1] = T(0);
//     coord[2] = T(0);
//     _coord[axis] = coord;
// }

template <typename T>
T clsVec<T>::x() const {return _coord[0];}
template <typename T>
T clsVec<T>::y() const {return _coord[1];}
template <typename T>
T clsVec<T>::z() const {return _coord[2];}

template <typename T>
T &clsVec<T>::x() {return _coord[0];}
template <typename T>
T &clsVec<T>::y() {return _coord[1];}
template <typename T>
T &clsVec<T>::z() {return _coord[2];}

template <typename T>
void clsVec<T>::x(T xi) {_coord[0] = xi;} 
template <typename T>
void clsVec<T>::y(T yi) {_coord[1] = yi;}
template <typename T>
void clsVec<T>::z(T zi) {_coord[2] = zi;}

template <typename T>
clsVec<T> clsVec<T>::operator+(const clsVec<T>& other) const {
    return {_coord[0] + other.x(), _coord[1] + other.y(), _coord[2] + other.z()};
}  

template <typename T>
clsVec<T> clsVec<T>::operator-(const clsVec<T>& other) const {
    return {_coord[0] - other.x(), _coord[1] - other.y(), _coord[2] - other.z()};
}  

template <typename T>
clsVec<T> clsVec<T>::operator+() const {
    return {_coord[0], _coord[1], _coord[2]};
}  

template <typename T>
clsVec<T> clsVec<T>::operator-() const {
    return {-_coord[0], -_coord[1], -_coord[2]};
}

template <typename T>
T& clsVec<T>::operator[](const int& idx){
    return _coord[idx];
}

template class clsVec<size_t>;
template class clsVec<float>;
template class clsVec<double>;
