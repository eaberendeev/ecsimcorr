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

//typedef Eigen::SparseMatrix<double, Eigen::ColMajor> Operator;

struct IndVal{
    int ind;
    double val;
};

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

    T& x() {return d[0];}
    T x() const {return d[0];}
    T& y() {return d[1];}
    T y() const {return d[1];}
    T& z() {return d[2];}
    T z() const {return d[2];}
    T& r() {return d[1];}
    T r() const {return d[1];}
    T& p() {return d[2];}
    T p() const {return d[2];}
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

//output
template<typename T>    //ostream output
std::ostream& operator<<(std::ostream &out, const vec3<T>& v) {
    out<<v(0)<<" "<<v(1)<<" "<<v(2);
    return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;
using int3 = vec3<int>;
template <typename T> 
struct Array3D{

    Array3D(int n1, int n2, int n3){
       allocate(n1,n2,n3);
    }
    Array3D(const int3& nn){
       allocate(nn(0),nn(1),nn(2));
    }
    Array3D(Array3D &&other): _data{other._data}, 
                  _size1{other._size1},_size2{other._size2},_size3{other._size3} {
        other._data = nullptr;
        other._size1 = other._size2 = other._size3 =  0; 
    }
    Array3D(const Array3D &other): _size1{other._size1},_size2{other._size2},_size3{other._size3} { 
       	allocate( _size1,_size2,_size3 );
    	for(auto i = 0; i < capacity(); ++i){
            _data[i] = other._data[i] ;
        }
    }
    Array3D(){
        _data = nullptr;
        _size1 = _size2 =_size3 = 0.;
    }
    
    void allocate(int n1, int n2, int n3){
        _data = new T[ n1 * n2 * n3];
        _size1 = n1;
        _size2 = n2;
        _size3 = n3;
    }
    
    void clear(){
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = 0.;
        }
    }

    void free(){
        if (_data != nullptr)
            delete[] _data;
        _size1 = _size2 =_size3 = 0.;
    }

    ~Array3D(){
       free();
    }
    
    Array3D<T>& operator=(T s) {
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = s ;
        }
        return (*this);
    }
    
    Array3D<T>& operator=(const Array3D<T>& array3D) {
    	// Self-assignment check
    	if (this == &array3D)
        	return *this;
        // Checking the conformity of dimensions
       	assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = array3D._data[i] ;
        }
        return (*this);
    }   
    Array3D<T>& operator=(Array3D<T> &&array3D) {
    	// Self-assignment check
    	if (this == &array3D)
        	return *this;
        // Checking the conformity of dimensions
       	assert( capacity() == array3D.capacity() );
		delete[] _data;
       	_data = array3D._data;
        _size1 = array3D._size1;
        _size2 = array3D._size2;
        _size3 = array3D._size3;
        array3D._data = nullptr;
        array3D._size1 = array3D._size2 = array3D._size3 =  0; 
        return (*this);
    }   

    Array3D<T>& operator-=(const Array3D<T>& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= array3D._data[i] ;
        }
        return (*this);
    }
    
    Array3D<T>& operator+=(const Array3D<T>& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += array3D._data[i] ;
        }
        return (*this);
    }          
    
    T& operator() (int i, int j, int k) {
        #if DEBUG > 0 
            if ( capacity() <= i * _size2 * _size3 + j * _size3 + k || i * _size2 * _size3 + j * _size3 + k < 0 ){
                std::string msg = "IndexError3";
                std::cout << msg << ": "<< int3(i,j,k)<< "\n";
                throw msg;
            };
        #endif
    	return _data[i * _size2 * _size3 + j * _size3 + k ];
    }

    const T& operator() (int i, int j, int k) const{ 
        #if DEBUG > 0 
            if ( capacity() <= i * _size2 * _size3 + j * _size3 + k || i * _size2 * _size3 + j * _size3 + k < 0 ){
                std::string msg = "IndexError3";
                std::cout << msg << ": "<< int3(i,j,k) << "\n";
                throw msg;
            };
        #endif
    	return _data[i * _size2 * _size3 + j * _size3 + k];
    }
    T& operator() (int i) {
    	return _data[i];
    }
    const T& operator() (int i) const {
    	return _data[i];
    }
    T& data(int i) {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0 ){
                std::string msg = "IndexError3";
                std::cout << msg << ": " << i << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
    const T& data(int i) const {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0){
                std::string msg = "IndexError3";
                std::cout << msg << ": " << i << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
  
    int3 size() const{
        return int3(_size1,_size2,_size3);
    }
    int capacity() const{
        return _size1*_size2*_size3;
    }

protected:
    T* _data;
    int _size1, _size2, _size3;
};

template <typename T> 
struct Array2D{
    Array2D(int n1, int n2){
       allocate(n1,n2);
    }
    Array2D(int2 nn){
       allocate(nn(0),nn(1));
    }
    Array2D(Array2D &&other): _data{other._data}, _size1{other._size1}, _size2{other._size2}{
        other._data = nullptr;
        other._size1 = other._size2 = 0; 
    }
    Array2D(){
        _data = nullptr;
        _size1 = _size2 = 0; 
    }
    
    void allocate(int n1, int n2){
        _data = new T[ n1 * n2];
        _size1 = n1;
        _size2 = n2;
    }
    
    void clear(){
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = 0.;
        }
    }
    void free(){
        if (_data != nullptr)
            delete[] _data;
        _size1 = _size2 = 0.;
    }
    
    ~Array2D(){
       free();
    }
    Array2D<T>& operator=(T s) {
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = s ;
        }
        return (*this);
    }

    Array2D<T>& operator=(const Array2D<T>& array2D) {
        assert( capacity() == array2D.capacity() && "Array 2D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = array2D._data[i] ;
        }
        return (*this);
    }
    Array2D<T>& operator=(const Array2D<T>&& other) {
        //assert( capacity() == other.capacity() );
        
        if(&other == this) return (*this);

        free();
        _size1 = other._size1;
        _size2 = other._size2;
        _data = other._data;
        other._data = nullptr;
        other._size1 = other._size2 = 0;
        
        return (*this);
    }

    Array2D<T>& operator-=(const Array2D<T>& array2D) {
        assert( capacity() == array2D.capacity()  && "Array 2D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= array2D._data[i] ;
        }
        return (*this);
    }
    
    Array2D<T>& operator+=(const Array2D<T>& array2D) {
        assert( capacity() == array2D.capacity()  && "Array 2D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += array2D._data[i] ;
        }
        return (*this);
    }

    T& operator() (int i, int j) {
        #if DEBUG > 0 
            if ( capacity() <= i * _size2  + j ||  i * _size2  + j < 0 ){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif
	return _data[i * _size2 + j];
    }

    const T& operator() (int i, int j) const{	
        #if DEBUG > 0 
            if ( capacity() <= i * _size2  + j ||  i * _size2  + j < 0 ){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif	
        return _data[i * _size2 + j];
    }
    
    T& data(int i) {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0 ){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
    const T& data(int i) const {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
  
    int2 size() const{
        return int2(_size1,_size2);
    }
    int capacity() const{
        return _size1*_size2;
    }
    T sum_d2(int i) const{
        T sum = 0;
        for( int t =  0; t < _size2; ++t){
          sum += _data[i * _size2 + t];
        }
        return sum;
    }

protected:
    T* _data;
    int _size1,_size2;
};

template <typename T> 
struct Array1D{
    
    Array1D(int size){
       allocate(size);
    }
    Array1D(Array1D &&other): _data{other._data}, _size{other._size} {
        other._data = nullptr;
        other._size = 0; 
    }
    Array1D(){
        _size = 0.;
    }
    
    void allocate(int sizeDim){
        _data = new T[ sizeDim];
        _size =  sizeDim;
    }
    
    void clear(){
        for(auto i = 0; i <_size; ++i){
            _data[i] = 0.;
        }
    }
    void free(){
        if (_data != nullptr)
            delete[] _data;
        _size = 0.;
    }
    
    Array1D<T>& operator=(T s) {
        for(auto i = 0; i < _size; ++i){
            _data[i] = s ;
        }
        return (*this);
    }

    Array1D<T>& operator=(const Array1D<T>& other) {
        assert( capacity() == other.capacity() && "Array 1D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = other._data[i] ;
        }
        return (*this);
    }

    Array1D<T>& operator-=(const Array1D<T>& other) {
        assert( capacity() == other.capacity()  && "Array 1D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= other._data[i] ;
        }
        return (*this);
    }
    
    Array1D<T>& operator+=(const Array1D<T>& other) {
        assert( capacity() == other.capacity()  && "Array 1D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += other._data[i] ;
        }
        return (*this);
    }
    Array1D<T>& operator=(const Array1D<T>&& other) {
        
        if( &other == this) return (*this);

        free();
        _size = other._size;
        _data = other._data;
        other._data = nullptr;
        other._size = 0;
        return (*this);
    }
    ~Array1D(){
        free();
    }
    
    T& operator() (int i) {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << ": index = " << i <<". capacity = "<< capacity() << "\n";
                throw msg;
            };
        #endif      
        return _data[i];
    }

    const T& operator() (int i) const{ 
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif          
        return _data[i];
    }
    int size() const{
        return _size;
    }
    int capacity() const{
        return _size;
    }
    T& data(int i) {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        return _data[i];
    }
    const T& data(int i) const {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        return _data[i];
    }  

protected:
    T *_data;
    int _size;  
};
/*
template <typename T> 
struct Array{

    Array(int maxSize){
        allocate(maxSize);
        _size = 0;
    }
    
    Array(){
        _data = nullptr;
        _maxSize = 0;
        _size = 0;
    };
    Array(Array &&arr): _data(arr._data), _size(arr._size), _maxSize(arr._maxSize) {
        arr._data = nullptr;
        arr._size = 0;
        arr._maxSize = 0;
    }
    ~Array(){
        free();
    }
        
    void allocate(int maxSize){
        _data = new T[maxSize];
        _maxSize = maxSize;
    }

    void clear(){
        _size = 0;
    }   
    
    void free(){
        if( _data!= nullptr)
            delete[] _data;
        _size = 0;
    }   
    

    void push_back(const T& elem){
        _data[_size] = elem;
        ++_size;
        #if DEBUG > 0 
            if ( _size >= capacity() ){
                std::string msg = "AddError";
                std::cout << msg << ": " << _size << " " << _maxSize << "\n";
                throw msg;
            };
        #endif  
    }
    
    void del(int k){
        #if DEBUG > 0 
            if ( _size < k ||  k < 0 ){
                std::string msg = "DelError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        --_size;
        _data[k] = _data[_size];
    }
    int size() const{
        return _size;
    }
    void resize(int newSize){
        _size = newSize;
    }
    int capacity() const{
        return _maxSize;
    }
    T& operator() (int i) {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexErrorP";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif          
        return _data[i];
    }

    const T& operator() (int i) const{
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexErrorP";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        return _data[i];
    }
    T back(){
        return _data[_size-1];
    }
    void pop_back(){
        --_size;
    }
        T* _data;
        int _size, _maxSize;
};
*/

struct Field3d{

    Field3d(int n1, int n2, int n3, int d){
       allocate(n1,n2,n3,d);
    }
    Field3d(int3 nn, int d){
       allocate(nn(0),nn(1),nn(2),d);
    }

    Field3d(Field3d &&other): _data{other._data}, 
                  _size1{other._size1},_size2{other._size2},_size3{other._size3}, _nd{other._nd} {
    }
    Field3d(const Field3d &other): _size1{other._size1},_size2{other._size2},_size3{other._size3}, _nd{other._nd} { 
        allocate( _size1,_size2,_size3,_nd );
            _data = other._data ;
    }
    Field3d(){
        //_data = nullptr;
        _size1 = _size2 =_size3 = 0.;
    }
    
    void allocate(int n1, int n2, int n3, int d){
        _data.resize(n1* n2* n3* d);
        _size1 = n1;
        _size2 = n2;
        _size3 = n3;
        _nd = d;
    }
    
    void clear(){
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = 0.;
        }
    }

    // void free(){
    //     if (_data != nullptr)
    //         delete[] _data;
    //     _size1 = _size2 =_size3 = 0.;
    // }

    // ~Field3d(){
    //    //free();
    // }
    
    Field3d& operator=(double s) {
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = s ;
        }
        return (*this);
    }
    
    Field3d& operator=(const Field3d& array3D) {
        // Self-assignment check
        if (this == &array3D)
            return *this;
        // Checking the conformity of dimensions
        assert( capacity() == array3D.capacity() );
        //for(auto i = 0; i < capacity(); ++i){
            _data = array3D._data ;
        //}
        return (*this);
    }   
    Field3d& operator=(Field3d &&array3D) {
        // Self-assignment check
        if (this == &array3D)
            return *this;
        // Checking the conformity of dimensions
        assert( capacity() == array3D.capacity() );
        //delete[] _data;
        _data = array3D._data;
        _size1 = array3D._size1;
        _size2 = array3D._size2;
        _size3 = array3D._size3;
        _nd = array3D._nd;
        //array3D._data = nullptr;
        //array3D._size1 = array3D._size2 = array3D._size3 =  0; 
        return (*this);
    }   

    Field3d& operator-=(const Field3d& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= array3D._data[i] ;
        }
        return (*this);
    }
    
    Field3d& operator+=(const Field3d& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += array3D._data[i] ;
        }
        return (*this);
    }         


    double& operator() (int i, int j, int k, int d) {
        return _data[d + _nd*(i * _size2 * _size3 + j * _size3 + k) ];
    }
    const double& operator() (int i, int j, int k, int d) const{
        return _data[d + _nd*(i * _size2 * _size3 + j * _size3 + k) ];
    }

    double& operator() (int i) {
        return _data[i];
    }
    const double& operator() (int i) const{
        return _data[i];
    }

    int nd() const{
        return _nd;
    }

    int3 size() const{
        return int3(_size1,_size2,_size3);
    }
    int capacity() const{
        return _nd*_size1*_size2*_size3;
    }

    Eigen::VectorXd& data() { 
        return _data;
    }
    const Eigen::VectorXd& data() const {
        return _data;
    }  

    Eigen::VectorXd _data;
    int _size1, _size2, _size3, _nd;
};


struct Field2d{

    Field2d(int n1, int n2,int d){
       allocate(n1,n2,d);
    }

    Field2d(Field2d &&other): _data{other._data}, 
                  _size1{other._size1},_size2{other._size2}, _nd{other._nd} {
    }
    Field2d(const Field2d &other): _size1{other._size1},_size2{other._size2}, _nd{other._nd} { 
        allocate( _size1,_size2,_nd );
            _data = other._data ;
    }
    Field2d(){
        //_data = nullptr;
        _size1 = _size2 = 0.;
    }
    
    void allocate(int n1, int n2, int d){
        _data.resize(n1* n2* d);
        _size1 = n1;
        _size2 = n2;
        _nd = d;
    }
    
    void clear(){
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = 0.;
        }
    }

    // void free(){
    //     if (_data != nullptr)
    //         delete[] _data;
    //     _size1 = _size2 =_size3 = 0.;
    // }

    // ~Field3d(){
    //    //free();
    // }
    
    Field2d& operator=(double s) {
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = s ;
        }
        return (*this);
    }
    
    Field2d& operator=(const Field2d& array3D) {
        // Self-assignment check
        if (this == &array3D)
            return *this;
        // Checking the conformity of dimensions
        assert( capacity() == array3D.capacity() );
        //for(auto i = 0; i < capacity(); ++i){
            _data = array3D._data ;
        //}
        return (*this);
    }   
    Field2d& operator=(Field2d &&array3D) {
        // Self-assignment check
        if (this == &array3D)
            return *this;
        // Checking the conformity of dimensions
        assert( capacity() == array3D.capacity() );
        //delete[] _data;
        _data = array3D._data;
        _size1 = array3D._size1;
        _size2 = array3D._size2;
        _nd = array3D._nd;
        //array3D._data = nullptr;
        //array3D._size1 = array3D._size2 = array3D._size3 =  0; 
        return (*this);
    }   

    Field2d& operator-=(const Field3d& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= array3D._data[i] ;
        }
        return (*this);
    }
    
    Field2d& operator+=(const Field2d& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += array3D._data[i] ;
        }
        return (*this);
    }         


    double& operator() (int i, int j, int d) {
        return _data[d + _nd*(i * _size2 + j) ];
    }
    const double& operator() (int i, int j, int d) const{
        return _data[d + _nd*(i * _size2  + j ) ];
    }

    int nd() const{
        return _nd;
    }

    int2 size() const{
        return int2(_size1,_size2);
    }
    int capacity() const{
        return _nd*_size1*_size2;
    }

    Eigen::VectorXd _data;
    int _size1, _size2, _nd;
};

template <class T>
struct Array : std::vector<T> {
    using std::vector<T>::reserve;
    //using std::vector<T>::size;
    using std::vector<T>::pop_back;
    using std::vector<T>::erase;
    Array(int maxSize){
        reserve(maxSize);
        //_size = 0;
    }    
    T& operator() (int i) {
        return (*this)[i];
    }

    const T& operator() (int i) const{
        return (*this)[i];
    }
    void del(int k){
        (*this)[k] = (*this)[ size()-1];
        pop_back();
    }
    int size() const{
        return int(std::vector<T>::size() );
    }

};
#endif 
