/*
    This file is part of LibImportance library that provides a technique for guiding
    transport paths towards the important places in the scene. This is a direct implementation
    of the method described in the paper "On-line Learning of Parametric Mixture 
    Models for Light Transport Simulation", ACM Trans. Graph. (SIGGRAPH 2014) 33, 4 (2014).
   
    Copyright (c) 2014 by Jiri Vorba, Ondrej Karlik, Martin Sik.

    LibImportance library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    LibImportance library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once

#include "Sse.h"
#include "Vector3.h"

namespace Importance {
    
template<class T>
class MultiVector3{
public:

    T x, y, z;
    

    IMPORTANCE_INLINE MultiVector3() { }
    IMPORTANCE_INLINE MultiVector3(const T& x, const T& y, const T& z) : x(x), y(y), z(z) { }
    IMPORTANCE_INLINE MultiVector3(const float x, const float y, const float z) : x(x), y(y), z(z) { }
    IMPORTANCE_INLINE explicit MultiVector3(const Vector3 vect) : x(vect.x), y(vect.y), z(vect.z) { }
    
    //IMPORTANCE_INLINE Vector3 operator[](const int i) const {
    //    return Vector3(x[i], y[i], z[i]);
    //}
    //} 
    IMPORTANCE_INLINE Vector3 get(const int i) const {
        return Vector3(x[i], y[i], z[i]);
    }
    IMPORTANCE_INLINE void set(const int i, const Vector3 input) {
        x[i] = input.x;
        y[i] = input.y;
        z[i] = input.z;
    }

    IMPORTANCE_INLINE MultiVector3 operator+(const MultiVector3& vect) const {
        return MultiVector3(vect.x+x, vect.y+y, vect.z+z);
    }
    IMPORTANCE_INLINE MultiVector3 operator-(const MultiVector3& vect) const {
        return MultiVector3(x-vect.x, y-vect.y, z-vect.z);
    }
    IMPORTANCE_INLINE MultiVector3 operator+(const Vector3 vect) const {
        return MultiVector3(vect.x()+x, vect.y()+y, vect.z()+z);
    }
    IMPORTANCE_INLINE MultiVector3 operator*(const T& vect) const {
        return MultiVector3(vect*x, vect*y, vect*z);
    }
    IMPORTANCE_INLINE MultiVector3& operator*=(const T& vect) {
        x *= vect;
        y *= vect;
        z *= vect;
        return *this;
    }
    IMPORTANCE_INLINE MultiVector3& operator/=(const T& vect) {
        const T inv = vect.getInverse();
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }

    IMPORTANCE_INLINE MultiVector3 operator*(const float vect) const {
        return MultiVector3(vect*x, vect*y, vect*z);
    }
    IMPORTANCE_INLINE MultiVector3& operator*=(const Vector3 vect) {
        x *= vect.x();
        y *= vect.y();
        z *= vect.z();
        return *this;
    }
    IMPORTANCE_INLINE MultiVector3& operator+=(const MultiVector3& vect) {
        x += vect.x;
        y += vect.y;
        z += vect.z;
        return *this;
    }
    IMPORTANCE_INLINE MultiVector3& operator+=(const Vector3 vect)  {
        this->x += vect.x();
        this->y += vect.y();
        this->z += vect.z();
        return *this;
    }
    IMPORTANCE_INLINE MultiVector3 operator*=(const float val) {
        x *= val;
        y *= val;
        z *= val;
        return *this;
    }
    IMPORTANCE_INLINE T squares() const {
        return x*x + y*y + z*z;
    }
    IMPORTANCE_INLINE T sizesFast() const {
        return squares().sqrtFast();
    }
    IMPORTANCE_INLINE T sizes() const {
        return squares().sqrt();
    }

    IMPORTANCE_INLINE Vector3 avg() const {
        return Vector3(x.avg(), y.avg(), z.avg());
    }
    IMPORTANCE_INLINE MultiVector3 abs() const {
        return MultiVector3(x.abs(), y.abs(), z.abs());
    }
    IMPORTANCE_INLINE Float4 max() const {
        return Float4::max(x, Float4::max(y, z));
    }
    IMPORTANCE_INLINE MultiVector3 clamp(const float minimum, const float maximum) const {
        return MultiVector3(x.clamp(minimum, maximum), y.clamp(minimum, maximum),z.clamp(minimum, maximum));
    }

    IMPORTANCE_INLINE friend MultiVector3<T> operator*(const float factor, const MultiVector3<T>& vect) {
        return MultiVector3<T>(factor*vect.x, factor*vect.y, factor*vect.z);
    }

    IMPORTANCE_INLINE static T dot(const MultiVector3<T>& first, const MultiVector3<T>& second) {
        return first.x*second.x + first.y*second.y + first.z*second.z;
    }
    IMPORTANCE_INLINE static T dot(const MultiVector3<T>& first, const Vector3 second) {
        return first.x*second.x + first.y*second.y + first.z*second.z;
    }
    IMPORTANCE_INLINE bool isNormalized() const {
        return  ((this->sizes() - Float4(1.f)).abs() < Float4(1e-4f)).allTrue();
    }
};

//IMPORTANCE_INLINE MultiVector3<T> operator*(const Vector3 vect, const T& factors) {
//    return MultiVector3<T>(vect.x()*factors, vect.y()*factors, vect.z()*factors);
//}
//template<class T>
//IMPORTANCE_INLINE MultiVector3<T> operator+(const Vector3 offset, const MultiVector3<T>& vect) {
//    return MultiVector3<T>(offset.x()+vect.x, offset.y()+vect.y, offset.z()+vect.z);
//}
//template<class T>
//IMPORTANCE_INLINE MultiVector3<T> operator*(const Vector3 mult, const MultiVector3<T>& vect) {
//    return MultiVector3<T>(mult.x()*vect.x, mult.y()*vect.y, mult.z()*vect.z);
//}

typedef MultiVector3<Float4> QuadVector3;
typedef MultiVector3<Float8> OctoVector3;



/*
class QuadVector3d : public Corona::Object {
public:

    Double4 x, y, z;


    IMPORTANCE_INLINE QuadVector3d() { }
    IMPORTANCE_INLINE QuadVector3d(const Double4& x, const Double4& y, const Double4& z) : x(x), y(y), z(z) { }
    IMPORTANCE_INLINE QuadVector3d(const double x, const double y, const double z) : x(y), y(y), z(z) { }
    IMPORTANCE_INLINE QuadVector3d operator+(const QuadVector3d& vect) const {
        return QuadVector3d(vect.x+x, vect.y+y, vect.z+z);
    }

    IMPORTANCE_INLINE QuadVector3d operator*(const Double4& vect) const {
        return QuadVector3d(vect*x, vect*y, vect*z);
    }

    IMPORTANCE_INLINE QuadVector3d operator+=(const QuadVector3d& vect) {
        x += vect.x;
        y += vect.y;
        z += vect.z;
        return *this;
    }

    IMPORTANCE_INLINE Double4 squares() const {
        return x*x + y*y + z*z;
    }
    IMPORTANCE_INLINE Double4 sizes() const {
        return squares().sqrt();
    }
};

*/



}