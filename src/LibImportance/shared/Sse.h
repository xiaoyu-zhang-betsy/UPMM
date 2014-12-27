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
/*<~API_TAG~>*/

#pragma warning(push)
#pragma warning(disable:4201)

#include "FastFloat/FastFloat.h"
#include <xmmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>
#include "FastFloat/avx_mathfun.h"

namespace Importance {

namespace Sse {
    template<int a, int b, int c, int d>
    IMPORTANCE_INLINE __m128 shuffle(const __m128 in) {
        return _mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(in), _MM_SHUFFLE(d, c, b, a)));
    }
    template<int a, int b, int c, int d>
    IMPORTANCE_INLINE __m128i shuffle(const __m128i in) {
        return _mm_shuffle_epi32(in, _MM_SHUFFLE(d, c, b, a));
    }

    IMPORTANCE_INLINE float get(const __m128 in, const int index) {
        IMPORTANCE_ASSERT(unsigned(index) < 4);
        return ((float*)(&in))[index];
    }
    IMPORTANCE_INLINE float& getRef(__m128& in, const int index) {
        IMPORTANCE_ASSERT(unsigned(index) < 4);
        return ((float*)(&in))[index];
    }
    IMPORTANCE_INLINE float get(const __m256& in, const int index) {
        IMPORTANCE_ASSERT(unsigned(index) < 8);
        return ((float*)(&in))[index];
    }
    IMPORTANCE_INLINE float& getRef(__m256& in, const int index) {
        IMPORTANCE_ASSERT(unsigned(index) < 8);
        return ((float*)(&in))[index];
    }


    //IMPORTANCE_INLINE double getFromAvx(const __m256d& in, const int index) {
    //    CASSERT(unsigned(index) < 4);
    //    return ((double*)(&in))[index];
    //}
    //IMPORTANCE_INLINE double& getRefFromAvx(__m256d& in, const int index) {
    //    CASSERT(unsigned(index) < 4);
    //    return ((double*)(&in))[index];
    //}
}


class Int4;
class Float4;

class Bool4 {
protected:

    __declspec(align(16)) union {
        __m128i _sse;
        int _data[4];
    };

    IMPORTANCE_INLINE bool isValid() const {
        for(int i = 0; i < 4; i++) {
            if((_data[i] != 0 && _data[i] != -1)) {
                return false;
            }
        }
        return true;
    } 

public:

    IMPORTANCE_INLINE explicit Bool4(const __m128i _sse) : _sse(_sse) { }
    IMPORTANCE_INLINE explicit Bool4(const __m128 _sse)  : _sse(_mm_castps_si128(_sse)) {
        IMPORTANCE_ASSERT(isValid());
    }
    IMPORTANCE_INLINE explicit Bool4(const bool value) : _sse(_mm_set1_epi32(-1*value)) {
        IMPORTANCE_ASSERT(isValid());
    }

    IMPORTANCE_INLINE const bool operator[](const int index) const {
        IMPORTANCE_ASSERT(unsigned(index) <= 3 && isValid());
        return this->_data[index] != 0;
    }


    IMPORTANCE_INLINE Bool4 operator&&(const Bool4 other) const {
        return Bool4(_mm_and_si128(_sse, other._sse));
    }
    IMPORTANCE_INLINE Bool4 operator||(const Bool4 other) const {
        return Bool4(_mm_or_si128(_sse, other._sse));
    }
    IMPORTANCE_INLINE Bool4 operator==(const Bool4 other) const {
        return Bool4(_mm_cmpeq_epi32(_sse, other._sse));
    }
    IMPORTANCE_INLINE bool allFalse() const {
        IMPORTANCE_ASSERT(isValid());
        return _mm_movemask_ps(_mm_castsi128_ps(_sse)) == 0;
    }
    IMPORTANCE_INLINE bool allTrue() const {
        IMPORTANCE_ASSERT(isValid());
        return _mm_movemask_ps(_mm_castsi128_ps(_sse)) == 0xF;
    }
    IMPORTANCE_INLINE bool anyFalse() const {
        return !allTrue();
    }
    IMPORTANCE_INLINE Int4 maskedInts(const unsigned int mask) const;
    
    IMPORTANCE_INLINE Float4 blend(const Float4 trueVals, const Float4 falseVals) const;
    IMPORTANCE_INLINE Int4 blend(const Int4 trueVals, const Int4 falseVals) const;

};


class Float4 {    
public:
    __m128 data;
    
    IMPORTANCE_INLINE Float4() {
#ifdef LIBIMP_DEBUG
        data = _mm_set1_ps(NAN);
#endif
    }
    IMPORTANCE_INLINE explicit Float4(const __m128 _sse)                                 : data(_sse) { }
    IMPORTANCE_INLINE Float4(const float x, const float y, const float z, const float w) : data(_mm_set_ps(w, z, y, x)) { }
    IMPORTANCE_INLINE explicit Float4(const float value)                                 : data(_mm_set1_ps(value)) { }
    IMPORTANCE_INLINE Float4(const float* memoryLocation)                                : data(_mm_load_ps(memoryLocation)) {
        IMPORTANCE_ASSERT(((__int64)memoryLocation%16) == 0);
    }

    IMPORTANCE_INLINE const float operator[](const int index) const {
       return Sse::get(data, index);
    }
    IMPORTANCE_INLINE float& operator[](const int index) {
       return Sse::getRef(data, index);
    }
    IMPORTANCE_INLINE Float4 operator+(const Float4 other) const {
        return Float4(_mm_add_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float4 operator-(const Float4 other) const {
        return Float4(_mm_sub_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float4 operator*(const Float4 other) const {
        return Float4(_mm_mul_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float4 operator/(const Float4 other) const {
        return Float4(_mm_div_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float4 operator+(const float fact) const {
        return Float4(_mm_add_ps(data, _mm_set1_ps(fact)));
    }
    IMPORTANCE_INLINE Float4 operator-(const float fact) const {
        return Float4(_mm_sub_ps(data, _mm_set1_ps(fact)));
    }
    IMPORTANCE_INLINE Float4 operator-() const {
        return Float4(_mm_sub_ps(_mm_set_ps1(0.f), data));
    }
    IMPORTANCE_INLINE Float4 operator*(const float fact) const {
        return Float4(_mm_mul_ps(data, _mm_set1_ps(fact)));
    }
    IMPORTANCE_INLINE Float4 operator/(const float fact) const {
        return Float4(_mm_div_ps(data, _mm_set1_ps(fact)));
    }
    friend IMPORTANCE_INLINE Float4 operator+(const float fact, const Float4 vect) {
        return Float4(_mm_add_ps(_mm_set1_ps(fact), vect.data));
    }
    friend IMPORTANCE_INLINE Float4 operator-(const float fact, const Float4 vect) {
        return Float4(_mm_sub_ps(_mm_set1_ps(fact), vect.data));
    }
    friend IMPORTANCE_INLINE Float4 operator*(const float fact, const Float4 vect) {
        return Float4(_mm_mul_ps(_mm_set1_ps(fact), vect.data));
    }

    IMPORTANCE_INLINE Float4& operator+=(const Float4 other) {
        this->data = _mm_add_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float4& operator-=(const Float4 other) {
        this->data = _mm_sub_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float4& operator*=(const Float4 other) {
        this->data = _mm_mul_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float4& operator/=(const Float4 other) {
        this->data = _mm_div_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float4& operator+=(const float other) {
        this->data = _mm_add_ps(data, _mm_set1_ps(other));
        return *this;
    }
    //IMPORTANCE_INLINE Float4 operator-() const {
    //    return Float4(_mm_xor_ps(data, _mm_castsi128_ps(_mm_set1_epi32(0x80000000))));
    //}
    IMPORTANCE_INLINE Float4 abs() const {
        return Float4(_mm_and_ps(data, _mm_castsi128_ps(_mm_set1_epi32(0x7fFFffFF))));
    }    
    IMPORTANCE_INLINE Float4& operator*=(const float factor) {
        this->data = _mm_mul_ps(data, _mm_set1_ps(factor));
        return *this;
    }
    IMPORTANCE_INLINE Float4& operator/=(const float factor) {
        this->data = _mm_div_ps(data, _mm_set1_ps(factor));
        return *this;
    }

    IMPORTANCE_INLINE Float4 getInverse() const {
        return Float4(_mm_div_ps(_mm_set1_ps(1.f), data));
    }
    IMPORTANCE_INLINE bool isReal() const {
        return (((*this) < INFINITY) && ((*this) > -INFINITY) && Bool4(_mm_cmpeq_ps(data, data))).allTrue();
    }

    static IMPORTANCE_INLINE Float4 min(const Float4 v1, const Float4 v2) {
        return Float4(_mm_min_ps(v1.data, v2.data));
    }
    static IMPORTANCE_INLINE Float4 max(const Float4 v1, const Float4 v2) {
        return Float4(_mm_max_ps(v1.data, v2.data));
    }

    IMPORTANCE_INLINE Bool4 operator<(const Float4 other) const {
        return Bool4(_mm_cmplt_ps(data, other.data));
    }
    IMPORTANCE_INLINE Bool4 operator>(const Float4 other) const {
        return Bool4(_mm_cmpgt_ps(data, other.data));
    }
    IMPORTANCE_INLINE Bool4 operator<=(const Float4 other) const {
        return Bool4(_mm_cmple_ps(data, other.data));
    }
    IMPORTANCE_INLINE Bool4 operator>=(const Float4 other) const {
        return Bool4(_mm_cmpge_ps(data, other.data));
    }
    IMPORTANCE_INLINE Bool4 operator<(const float other) const {
        return Bool4(_mm_cmplt_ps(data, _mm_set1_ps(other)));
    }
    IMPORTANCE_INLINE Bool4 operator>(const float other) const {
        return Bool4(_mm_cmpgt_ps(data, _mm_set1_ps(other)));
    }
    IMPORTANCE_INLINE Bool4 operator<=(const float other) const {
        return Bool4(_mm_cmple_ps(data, _mm_set1_ps(other)));
    }
    IMPORTANCE_INLINE Bool4 operator>=(const float other) const {
        return Bool4(_mm_cmpge_ps(data, _mm_set1_ps(other)));
    }
    IMPORTANCE_INLINE Bool4 operator!=(const Float4 other) const {
        return Bool4(_mm_cmpneq_ps(data, other.data));
    }
    

    IMPORTANCE_INLINE Float4 fastPow(const float exponent) const {
        return Float4(Ff::powFast(data, _mm_set1_ps(exponent)));
    }
    IMPORTANCE_INLINE Float4 log() const {
        return Float4(Ff::log(data));
    }
    IMPORTANCE_INLINE Float4 exp() const {
        return Float4(Ff::exp(data));
    }

    IMPORTANCE_INLINE Float4 sqrtFast() const {
        return Float4(_mm_rcp_ps(_mm_rsqrt_ps(data)));
    }
    IMPORTANCE_INLINE Float4 sqrt() const {
        return Float4(_mm_sqrt_ps(data));
    }
    IMPORTANCE_INLINE Float4 invSqrt() const {
        return sqrt().getInverse();
    }
    IMPORTANCE_INLINE Float4 invSqrtFast() const {
        return Float4(_mm_rsqrt_ps(data));
    }

    IMPORTANCE_INLINE float avg() const {  //TODO neni na to nejaky sse?
        return ((Sse::get(data, 0) + Sse::get(data, 1)) + (Sse::get(data, 2) + Sse::get(data, 3)))*0.25f;
    }    

    static IMPORTANCE_INLINE float dot(const Float4& v1, const Float4& v2) {
        return _mm_cvtss_f32(_mm_dp_ps(v1.data, v2.data, 0xF1));
    }
    IMPORTANCE_INLINE float square() const {
        return dot(*this, *this);
    }

	IMPORTANCE_INLINE float l2normFast() const {
		return Ff::sqrtFast(dot(*this, *this));
	}
    IMPORTANCE_INLINE float l2norm() const {
        return ::sqrt(dot(*this, *this));
    }
    IMPORTANCE_INLINE Float4 clamp(const float minimum, const float maximum) const {
        return Float4(_mm_min_ps(_mm_max_ps(_mm_set1_ps(minimum), data), _mm_set1_ps(maximum)));
    }

    static IMPORTANCE_INLINE Bool4 epsEqual(const Float4 x, const Float4 y) {
        return Bool4((x-y).abs() < 1e-6f);
    }
};

IMPORTANCE_INLINE std::ostream& operator<<(std::ostream& os, const Float4 sse) {
    os << sse[0] << ' ' << sse[1] << ' ' << sse[2] << ' ' << sse[3];
    return os;
}
IMPORTANCE_INLINE std::istream& operator>>(std::istream& is, Float4& sse) {
    is >> sse[0] >> sse[1] >> sse[2] >> sse[3];
    return is;
}

IMPORTANCE_INLINE bool isReal(const Float4 in) {
    return in.isReal();
}


class Int4 {
public:
    __declspec(align(16)) union {
        __m128i _sse;

        struct { 
            int x;
            int y;
            int z;
            int w;
        };

        int _data[4];
    };

    IMPORTANCE_INLINE Int4() { };
    IMPORTANCE_INLINE explicit Int4(const __m128i _sse) : _sse(_sse) { }
    IMPORTANCE_INLINE explicit Int4(const int value) : _sse(_mm_set1_epi32(value)) { }
    IMPORTANCE_INLINE Int4(const int x, const int y, const int z, const int w) : _sse(_mm_set_epi32(w, z, y, x)) { }
    IMPORTANCE_INLINE explicit Int4(const int* data) : _sse(_mm_castps_si128(_mm_load_ps((const float*)data))) {
        IMPORTANCE_ASSERT(((__int64)data%16) == 0);
    }

    static IMPORTANCE_INLINE Int4 reinterpretFloat(const Float4 input) {
        return Int4(_mm_castps_si128(input.data));
    }

    IMPORTANCE_INLINE const int operator[](const int index) const {
        IMPORTANCE_ASSERT(unsigned(index) <= 3);
        return this->_data[index];
    }
    IMPORTANCE_INLINE int& operator[](const int index) {
        IMPORTANCE_ASSERT(unsigned(index) <= 3);
        return this->_data[index];
    }

    IMPORTANCE_INLINE Int4 operator+(const Int4 other) const {
        return Int4(_mm_add_epi32(_sse, other._sse));
    }
    IMPORTANCE_INLINE Int4 operator-(const Int4 other) const {
        return Int4(_mm_sub_epi32(_sse, other._sse));
    }
    IMPORTANCE_INLINE Int4 operator*(const Int4 other) const {
        return Int4(_mm_mul_epi32(_sse, other._sse));
    }
    IMPORTANCE_INLINE Int4 operator+(const int fact) const {
        return Int4(_mm_add_epi32(_sse, _mm_set1_epi32(fact)));
    }
    IMPORTANCE_INLINE Int4 operator-(const int fact) const {
        return Int4(_mm_sub_epi32(_sse, _mm_set1_epi32(fact)));
    }
    IMPORTANCE_INLINE Int4 operator*(const int fact) const {
        return Int4(_mm_mul_epi32(_sse, _mm_set1_epi32(fact)));
    }
    IMPORTANCE_INLINE Int4& operator+=(const Int4 other) {
        this->_sse = _mm_add_epi32(_sse, other._sse);
        return *this;
    }
    IMPORTANCE_INLINE Int4& operator-=(const Int4 other) {
        _sse = _mm_sub_epi32(_sse, other._sse);
        return *this;
    }
    static IMPORTANCE_INLINE Int4 min(const Int4 x, const Int4 y) {
        return Int4(_mm_min_epi32(x._sse, y._sse));
    }

    IMPORTANCE_INLINE int min() const {
        const __m128i h = _mm_min_epi32(Sse::shuffle<1,0,3,2>(_sse), _sse);
        const int res = _mm_extract_epi32(_mm_min_epi32(Sse::shuffle<2,3,0,1>(h), h), 1);
        //CASSERT(res <= (*this)[0] && res <= (*this)[1] && res <= (*this)[2] && res <= (*this)[3]);
        return res;
    }


    IMPORTANCE_INLINE Float4 toFloat() const {
        return Float4(_mm_cvtepi32_ps(_sse));
    }

    IMPORTANCE_INLINE static Int4 float2int(const Float4 data) {
        return Int4(_mm_cvtps_epi32(data.data));
    }
};

//vraci 0 pro false, masked hodnotu pro 1
IMPORTANCE_INLINE Int4 Bool4::maskedInts(const unsigned int mask) const {
    return Int4(_mm_and_si128(_sse, _mm_set1_epi32(mask)));
}

IMPORTANCE_INLINE Float4 Bool4::blend(const Float4 trueVals, const Float4 falseVals) const {
    return Float4(_mm_blendv_ps(falseVals.data, trueVals.data, _mm_castsi128_ps(this->_sse)));
}
IMPORTANCE_INLINE Int4 Bool4::blend(const Int4 trueVals, const Int4 falseVals) const {
    return Int4(_mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(falseVals._sse), _mm_castsi128_ps(trueVals._sse), _mm_castsi128_ps(this->_sse))));
}



//AVX stuff



class Float8 {    
public:

    static const int DIMENSION = 8;

    __declspec(align(32)) __m256 data;
    
    IMPORTANCE_INLINE Float8() {
#ifdef LIBIMP_DEBUG
        data = _mm256_set1_ps(NAN);
#endif
    }
    IMPORTANCE_INLINE explicit Float8(const __m256& data)                                 : data(data) { }
    //IMPORTANCE_INLINE Float8(const float x, const float y, const float z, const float w) : data(_mm256_set_ps(w, z, y, x)) { }
    IMPORTANCE_INLINE explicit Float8(const float value)                                 : data(_mm256_set1_ps(value)) { }
    IMPORTANCE_INLINE Float8(const float* memoryLocation)                                : data(_mm256_load_ps(memoryLocation)) {
        IMPORTANCE_ASSERT(((__int64)memoryLocation%32) == 0);
    }

    IMPORTANCE_INLINE const float operator[](const int index) const {
        return Sse::get(data, index);
    }
    IMPORTANCE_INLINE float& operator[](const int index) {
        return Sse::getRef(data, index);
    }
    IMPORTANCE_INLINE Float8 operator+(const Float8& other) const {
        return Float8(_mm256_add_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float8 operator-(const Float8& other) const {
        return Float8(_mm256_sub_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float8 operator*(const Float8& other) const {
        return Float8(_mm256_mul_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float8 operator/(const Float8& other) const {
        return Float8(_mm256_div_ps(data, other.data));
    }
    IMPORTANCE_INLINE Float8 operator+(const float fact) const {
        return Float8(_mm256_add_ps(data, _mm256_set1_ps(fact)));
    }
    IMPORTANCE_INLINE Float8 operator-(const float fact) const {
        return Float8(_mm256_sub_ps(data, _mm256_set1_ps(fact)));
    }
    IMPORTANCE_INLINE Float8 operator*(const float fact) const {
        return Float8(_mm256_mul_ps(data, _mm256_set1_ps(fact)));
    }
    IMPORTANCE_INLINE Float8 operator/(const float fact) const {
        return Float8(_mm256_div_ps(data, _mm256_set1_ps(fact)));
    }
    friend IMPORTANCE_INLINE Float8 operator+(const float fact, const Float8& vect) {
        return Float8(_mm256_add_ps(_mm256_set1_ps(fact), vect.data));
    }
    friend IMPORTANCE_INLINE Float8 operator-(const float fact, const Float8& vect) {
        return Float8(_mm256_sub_ps(_mm256_set1_ps(fact), vect.data));
    }
    friend IMPORTANCE_INLINE Float8 operator*(const float fact, const Float8& vect) {
        return Float8(_mm256_mul_ps(_mm256_set1_ps(fact), vect.data));
    }

    IMPORTANCE_INLINE Float8& operator+=(const Float8& other) {
        this->data = _mm256_add_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float8& operator-=(const Float8& other) {
        this->data = _mm256_sub_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float8& operator*=(const Float8& other) {
        this->data = _mm256_mul_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float8& operator/=(const Float8& other) {
        this->data = _mm256_div_ps(data, other.data);
        return *this;
    }
    IMPORTANCE_INLINE Float8& operator+=(const float other) {
        this->data = _mm256_add_ps(data, _mm256_set1_ps(other));
        return *this;
    }
    IMPORTANCE_INLINE Float8 operator-() const {
        return Float8(_mm256_xor_ps(data, _mm256_castsi256_ps(_mm256_set1_epi32(0x80000000))));
    }
    IMPORTANCE_INLINE Float8 abs() const {
        return Float8(_mm256_and_ps(data, _mm256_castsi256_ps(_mm256_set1_epi32(0x7fFFffFF))));
    }    
    IMPORTANCE_INLINE Float8& operator*=(const float factor) {
        this->data = _mm256_mul_ps(data, _mm256_set1_ps(factor));
        return *this;
    }
    IMPORTANCE_INLINE Float8& operator/=(const float factor) {
        this->data = _mm256_div_ps(data, _mm256_set1_ps(factor));
        return *this;
    }

    IMPORTANCE_INLINE Float8 getInverse() const {
        return Float8(_mm256_div_ps(_mm256_set1_ps(1.f), data));
    }
    IMPORTANCE_INLINE bool isReal() const {
        return (_mm256_movemask_ps(_mm256_cmp_ps(data, _mm256_set1_ps(INFINITY), _CMP_LT_OQ)) &
            _mm256_movemask_ps(_mm256_cmp_ps(data, _mm256_set1_ps(-INFINITY), _CMP_GT_OQ)) &
            _mm256_movemask_ps(_mm256_cmp_ps(data, data, _CMP_EQ_UQ))) == 255;
    }

    static IMPORTANCE_INLINE Float8 min(const Float8& v1, const Float8& v2) {
        return Float8(_mm256_min_ps(v1.data, v2.data));
    }
    static IMPORTANCE_INLINE Float8 max(const Float8& v1, const Float8& v2) {
        return Float8(_mm256_max_ps(v1.data, v2.data));
    }

    IMPORTANCE_INLINE Float8 sqrtFast() const {
        return Float8(_mm256_rcp_ps(_mm256_rsqrt_ps(data)));
    }
    IMPORTANCE_INLINE Float8 sqrt() const {
        return Float8(_mm256_sqrt_ps(data));
    }
    IMPORTANCE_INLINE Float8 exp() const {
        return Float8(exp256_ps(data));
    }
    IMPORTANCE_INLINE Float8 invSqrt() const {
        return sqrt().getInverse();
    }
    IMPORTANCE_INLINE Float8 invSqrtFast() const {
        return Float8(_mm256_rsqrt_ps(data));
    }
    IMPORTANCE_INLINE Float8 clamp(const float minimum, const float maximum) const {
        return Float8(_mm256_min_ps(_mm256_max_ps(_mm256_set1_ps(minimum), data), _mm256_set1_ps(maximum)));
    }  
    //static IMPORTANCE_INLINE float dot(const Float8& v1, const Float8& v2) {
    //    return (float&)(_mm256_dp_ps(v1.data, v2.data, 0xF1));
    //}
};

}