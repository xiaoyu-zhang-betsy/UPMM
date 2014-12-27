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

#include <xmmintrin.h>
#include <smmintrin.h>
#include <xutility>
#include <sstream>
#include <string>
#include <iostream>
#include "NumberUtils.h"
#include "FastFloat/FastFloat.h"
#include "../shared/Config.h"

#undef min
#undef max


#pragma warning(push)
#pragma warning(disable:4201)   //nonstandard extension used : nameless struct/union

namespace Importance {

    template<int a, int b, int c, int d>
    IMPORTANCE_INLINE __m128 shuffle(const __m128 in) {
        return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(in), _MM_SHUFFLE(d, c, b, a)));
    }
    template<int a, int b, int c, int d>
    IMPORTANCE_INLINE __m128i shuffle(const __m128i in) {
        return _mm_shuffle_epi32(in, _MM_SHUFFLE(d, c, b, a));
    }



    class Vector3f {
    public:
        union {
            struct {
                float x;
                float y;
                float z;
                float extraData;
            };
            __m128 data;
        };

        IMPORTANCE_INLINE Vector3f() {
#ifdef LIBIMP_DEBUG
            data = _mm_set1_ps(FLOAT_NAN);
#endif
        }

        IMPORTANCE_INLINE void serialize( std::ostream & output ) const {
            output.write( (char *) &x, sizeof( float ) );
            output.write( (char *) &y, sizeof( float ) );
            output.write( (char *) &z, sizeof( float ) );
        }

        IMPORTANCE_INLINE void deserialize( std::istream & input ) {
            input.read( (char *) &x, sizeof( float ) );
            input.read( (char *) &y, sizeof( float ) );
            input.read( (char *) &z, sizeof( float ) );
        }

        IMPORTANCE_INLINE explicit Vector3f(const __m128& _sse) : data(_sse) { }    

        IMPORTANCE_INLINE Vector3f(const float x, const float y, const float z) {
            this->x = x;
            this->y = y;
            this->z = z;
        }

        IMPORTANCE_INLINE explicit Vector3f(const float x) : data(_mm_set1_ps(x)) { }

        IMPORTANCE_INLINE float operator[](const int index) const {
            IMPORTANCE_ASSERT(unsigned(index) < 3);
            return ((float*)&data)[index];
        }
        IMPORTANCE_INLINE float& operator[](const int index) {
            IMPORTANCE_ASSERT(unsigned(index) < 3);
            return ((float*)&data)[index];
        }

        IMPORTANCE_INLINE float size() const {
            return sqrt(square());
        }
        IMPORTANCE_INLINE float sizeApprox() const {
            return Ff::sqrtFast(square());
        }
        IMPORTANCE_INLINE Vector3f getNormalizedApprox(float& outSize) {
            outSize = sizeApprox();
            return *this / outSize;
        }

        IMPORTANCE_INLINE float length() const {
            return size();
        }

        IMPORTANCE_INLINE float square() const {
            return dot(*this, *this);
        }

        IMPORTANCE_INLINE Vector3f getNormalized() const {
            return *this / size();
        }

        IMPORTANCE_INLINE bool isNormalized() const {
            return ::abs(this->square() - 1.f) < 1e-4f;
        }    
        IMPORTANCE_INLINE Vector3f operator+(const Vector3f other) const {
            return Vector3f(_mm_add_ps(data, other.data));
        }
        IMPORTANCE_INLINE Vector3f operator-(const Vector3f other) const {
            return Vector3f(_mm_sub_ps(data, other.data));
        }
        IMPORTANCE_INLINE Vector3f operator*(const Vector3f other) const {
            return Vector3f(_mm_mul_ps(data, other.data));
        }
        IMPORTANCE_INLINE Vector3f operator/(const Vector3f other) const {
            return Vector3f(_mm_div_ps(data, other.data));
        }
        IMPORTANCE_INLINE Vector3f& operator+=(const Vector3f other) {
            this->data = _mm_add_ps(data, other.data);
            return *this;
        }
        IMPORTANCE_INLINE Vector3f& operator-=(const Vector3f other) {
            this->data = _mm_sub_ps(data, other.data);
            return *this;
        }
        IMPORTANCE_INLINE Vector3f& operator*=(const Vector3f other) {
            this->data = _mm_mul_ps(data, other.data);
            return *this;
        }
        IMPORTANCE_INLINE Vector3f& operator/=(const Vector3f other) {
            this->data = _mm_div_ps(data, other.data);
            return *this;
        }
        IMPORTANCE_INLINE Vector3f operator-() const {
            return Vector3f(_mm_xor_ps(data, _mm_castsi128_ps(_mm_set1_epi32(0x80000000))));
        }
        IMPORTANCE_INLINE Vector3f abs() const {
            return Vector3f(_mm_and_ps(data, _mm_castsi128_ps(_mm_set1_epi32(0x7fFFffFF))));
        }
        IMPORTANCE_INLINE Vector3f operator*(const float factor) const {
            return Vector3f(_mm_mul_ps(data, _mm_set1_ps(factor)));
        }
        IMPORTANCE_INLINE Vector3f operator/(const float factor) const {
            return Vector3f(_mm_div_ps(data, _mm_set1_ps(factor)));
        }
        IMPORTANCE_INLINE Vector3f& operator*=(const float factor) {
            this->data = _mm_mul_ps(data, _mm_set1_ps(factor));
            return *this;
        }
        IMPORTANCE_INLINE Vector3f& operator/=(const float factor) {
            this->data = _mm_div_ps(data, _mm_set1_ps(factor));
            return *this;
        }
        
        friend IMPORTANCE_INLINE Vector3f operator*(const float factor, const Vector3f vector);

        friend IMPORTANCE_INLINE Vector3f operator/(const float factor, const Vector3f vector) {
            return Vector3f(_mm_div_ps(_mm_set1_ps(factor), vector.data));
        }
        IMPORTANCE_INLINE float min() const {
            return std::min(std::min(x, y), z);
        }
        IMPORTANCE_INLINE float max() const {
            return std::max(std::max(x, y), z);
        }    
        IMPORTANCE_INLINE int argMin() const {
            return (x < y) ? ((x < z) ? 0 : 2) : ((y < z) ? 1 : 2);
        }
        IMPORTANCE_INLINE int argMax() const {
            return (x > y) ? ((x > z) ? 0 : 2) : ((y > z) ? 1 : 2);
        }
        IMPORTANCE_INLINE float l1Norm() const {
            const Vector3f absVals = this->abs();
            return absVals.x + absVals.y + absVals.z;
        }    
        IMPORTANCE_INLINE Vector3f getInverse() const {
            return Vector3f(_mm_div_ps(_mm_set1_ps(1.f), data));
        }

        IMPORTANCE_INLINE bool isReal() const {
            return (_mm_movemask_ps(_mm_and_ps(_mm_and_ps(
                _mm_cmplt_ps(data, _mm_set1_ps(FLOAT_INFINITY)),
                _mm_cmpgt_ps(data, _mm_set1_ps(-FLOAT_INFINITY))            
                ),
                _mm_cmpeq_ps(data, data))) & 0x7) == 0x7;
        }

        IMPORTANCE_INLINE bool operator==(const Vector3f other) const {
            return (_mm_movemask_ps(_mm_cmpeq_ps(data, other.data)) & 0x7) == 0x7;
        }

        friend IMPORTANCE_INLINE Vector3f cross(const Vector3f v1, const Vector3f v2) {
            return Vector3f(_mm_sub_ps(_mm_mul_ps(shuffle<1,2,0,3>(v1.data), shuffle<2,0,1,3>(v2.data)), 
                _mm_mul_ps(shuffle<2,0,1,3>(v1.data), shuffle<1,2,0,3>(v2.data))));
        }

        friend IMPORTANCE_INLINE float dot(const Vector3f v1, const Vector3f v2);

        static IMPORTANCE_INLINE Vector3f min(const Vector3f v1, const Vector3f v2) {
            return Vector3f(_mm_min_ps(v1.data, v2.data));
        }
        static IMPORTANCE_INLINE Vector3f max(const Vector3f v1, const Vector3f v2) {
            return Vector3f(_mm_max_ps(v1.data, v2.data));
        }

        IMPORTANCE_INLINE Vector3f clamp(const Vector3f low, const Vector3f high) const {
            IMPORTANCE_ASSERT(low.x <= high.x && low.y <= high.y && low.z < high.z);
            return Vector3f::min(Vector3f::max(*this, low), high);
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "[" << x << ", " << y << ", " << z << "]";
            return oss.str();
        }

        IMPORTANCE_INLINE Vector3f singlePrecision() const {
            return *this;
        }
    };

    IMPORTANCE_INLINE float dot(const Vector3f v1, const Vector3f v2) {
        return _mm_cvtss_f32(_mm_dp_ps(v1.data, v2.data, 0x71));
    }

    IMPORTANCE_INLINE Vector3f operator*(const float factor, const Vector3f vector) {
        return Vector3f(_mm_mul_ps(vector.data, _mm_set1_ps(factor)));
    }


    class Vector3d {
    public:

        union {
            struct {
                double x;
                double y;
                double z;
            };

            double data[3];
        };

        IMPORTANCE_INLINE void serialize( std::ostream & output ) const {
            output.write( (char *) &x, sizeof( double ) );
            output.write( (char *) &y, sizeof( double ) );
            output.write( (char *) &z, sizeof( double ) );
        }

        IMPORTANCE_INLINE void deserialize( std::istream & input ) {
            input.read( (char *) &x, sizeof( double ) );
            input.read( (char *) &y, sizeof( double ) );
            input.read( (char *) &z, sizeof( double ) );
        }

        IMPORTANCE_INLINE Vector3d() {
#ifdef LIBIMP_DEBUG
            x = y = z = DOUBLE_NAN;
#endif
        }

        IMPORTANCE_INLINE explicit Vector3d(const Vector3f in) : x(in.x), y(in.y), z(in.z) { }

        IMPORTANCE_INLINE Vector3d(const double x, const double y, const double z) {
            this->x = x;
            this->y = y;
            this->z = z;
        }

        IMPORTANCE_INLINE explicit Vector3d(const double x) : x(x), y(x), z(x) { }

        IMPORTANCE_INLINE Vector3f singlePrecision() const {
            return Vector3f(float(x), float(y), float(z));
        }

        IMPORTANCE_INLINE double operator[](const int index) const {
            IMPORTANCE_ASSERT(unsigned(index) < 3);
            return data[index];
        }
        IMPORTANCE_INLINE double& operator[](const int index) {
            IMPORTANCE_ASSERT(unsigned(index) < 3);
            return data[index];
        }

        IMPORTANCE_INLINE double size() const {
            return sqrt(square());
        }

        IMPORTANCE_INLINE double length() const {
            return size();
        }

        IMPORTANCE_INLINE double square() const {
            return dot(*this, *this);
        }

        IMPORTANCE_INLINE Vector3d getNormalized() const {
            return *this / size();
        }

        IMPORTANCE_INLINE bool isNormalized() const {
            return ::abs(this->square() - 1.) < 1e-4;
        }    
        IMPORTANCE_INLINE Vector3d operator+(const Vector3d other) const {
            return Vector3d(x+other.x, y+other.y, z+other.z);
        }
        IMPORTANCE_INLINE Vector3d operator-(const Vector3d other) const {
            return Vector3d(x-other.x, y-other.y, z-other.z);
        }
        IMPORTANCE_INLINE Vector3d operator*(const Vector3d other) const {
            return Vector3d(x*other.x, y*other.y, z*other.z);
        }
        IMPORTANCE_INLINE Vector3d operator/(const Vector3d other) const {
            return Vector3d(x/other.x, y/other.y, z/other.z);
        }
        IMPORTANCE_INLINE Vector3d& operator+=(const Vector3d other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }
        IMPORTANCE_INLINE Vector3d& operator-=(const Vector3d other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }
        IMPORTANCE_INLINE Vector3d& operator*=(const Vector3d other) {
            x *= other.x;
            y *= other.y;
            z *= other.z;
            return *this;
        }
        IMPORTANCE_INLINE Vector3d& operator/=(const Vector3d other) {
            x /= other.x;
            y /= other.y;
            z /= other.z;
            return *this;
        }
        IMPORTANCE_INLINE Vector3d operator-() const {
            return Vector3d(-x, -y, -z);
        }
        IMPORTANCE_INLINE Vector3d abs() const {
            return Vector3d(::abs(x), ::abs(y), ::abs(z));
        }
        IMPORTANCE_INLINE Vector3d operator*(const double factor) const {
            return Vector3d(x*factor, y*factor, z*factor);
        }
        IMPORTANCE_INLINE Vector3d operator/(const double factor) const {
            const double inv = 1./factor;
            return Vector3d(x*inv, y*inv, z*inv);
        }
        IMPORTANCE_INLINE Vector3d& operator*=(const double factor) {
            x *= factor;
            y *= factor;
            z *= factor;
            return *this;
        }
        IMPORTANCE_INLINE Vector3d& operator/=(const double factor) {
            const double inv = 1./factor;
            x *= inv;
            y *= inv;
            z *= inv;
            return *this;
        }
        
        friend IMPORTANCE_INLINE Vector3d operator*(const double factor, const Vector3d vector);

        friend IMPORTANCE_INLINE Vector3d operator/(const double factor, const Vector3d vector) {
            return Vector3d(factor/vector.x, factor/vector.y, factor/vector.z);
        }
        IMPORTANCE_INLINE double min() const {
            return std::min(std::min(x, y), z);
        }
        IMPORTANCE_INLINE double max() const {
            return std::max(std::max(x, y), z);
        }    
        IMPORTANCE_INLINE int argMin() const {
            return (x < y) ? ((x < z) ? 0 : 2) : ((y < z) ? 1 : 2);
        }
        IMPORTANCE_INLINE int argMax() const {
            return (x > y) ? ((x > z) ? 0 : 2) : ((y > z) ? 1 : 2);
        }
        IMPORTANCE_INLINE double l1Norm() const {
            const Vector3d absVals = this->abs();
            return absVals.x + absVals.y + absVals.z;
        }    
        IMPORTANCE_INLINE Vector3d getInverse() const {
            return 1. / *this;
        }

        IMPORTANCE_INLINE bool isReal() const {
            return x > -DOUBLE_INFINITY && y > -DOUBLE_INFINITY && z > -DOUBLE_INFINITY &&
                x < DOUBLE_INFINITY && y < DOUBLE_INFINITY && z < DOUBLE_INFINITY &&
                !isNaN(x) && !isNaN(y) && !isNaN(z);
        }

        IMPORTANCE_INLINE bool operator==(const Vector3d other) const {
            return x == other.x && y == other.y && z == other.z;
        }

        friend IMPORTANCE_INLINE Vector3d cross(const Vector3d v1, const Vector3d v2) {
            return Vector3d((v1.y * v2.z) - (v1.z * v2.y), 
                (v1.z * v2.x) - (v1.x * v2.z),
                (v1.x * v2.y) - (v1.y * v2.x));
        }

        friend IMPORTANCE_INLINE double dot(const Vector3d v1, const Vector3d v2);

        static IMPORTANCE_INLINE Vector3d min(Vector3d v1, const Vector3d v2) {
            if(v2.x < v1.x) {
                v1.x = v2.x;
            }
            if(v2.y < v1.y) {
                v1.y = v2.y;
            }
            if(v2.z < v1.z) {
                v1.z = v2.z;
            }
            return v1;
        }
        static IMPORTANCE_INLINE Vector3d max(Vector3d v1, const Vector3d v2) {
            if(v2.x > v1.x) {
                v1.x = v2.x;
            }
            if(v2.y > v1.y) {
                v1.y = v2.y;
            }
            if(v2.z > v1.z) {
                v1.z = v2.z;
            }
            return v1;
        }

        IMPORTANCE_INLINE Vector3d clamp(const Vector3d low, const Vector3d high) const {
            IMPORTANCE_ASSERT(low.x <= high.x && low.y <= high.y && low.z < high.z);
            return Vector3d::min(Vector3d::max(*this, low), high);
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "[" << x << ", " << y << ", " << z << "]";
            return oss.str();
        }
    };

    IMPORTANCE_INLINE double dot(const Vector3d v1, const Vector3d v2) {
        return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
    }

    IMPORTANCE_INLINE Vector3d operator*(const double factor, const Vector3d vector) {
        return Vector3d(vector.x*factor, vector.y*factor, vector.z*factor);
    }

}

#pragma warning(pop)