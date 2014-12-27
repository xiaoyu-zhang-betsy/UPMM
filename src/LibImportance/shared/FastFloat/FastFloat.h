#ifndef ___IMPORTANCE_FAST_FLOAT___
#define ___IMPORTANCE_FAST_FLOAT___


#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <math.h>


namespace Importance {

#define SSE_LIBRARY_BUILD
#define FF_INLINE __forceinline

// http://gruntthepeon.free.fr/ssemath/
namespace MathFun {
#include "sse_mathfun.h"
}

// http://code.google.com/p/ut-sse/source/browse/trunk/examples/particleFilter/sse/?r=2
namespace SseWrapper {
#include "sseMath.h"
}

// http://jrfonseca.blogspot.com/2008/09/fast-sse2-pow-tables-or-polynomials.html
namespace FastPow {
#include "pow.h"
#include <intrin.h>
}

#pragma warning(push)
#pragma warning(disable:4238)
#pragma warning(disable:4244)

#include "fastonebigheader.h"


    namespace Ff {

        FF_INLINE __m128 sqrtPrecise(const __m128 in) {
            return _mm_sqrt_ps(in);
        }
        FF_INLINE __m128 sqrtFast(const __m128 in) {
            return _mm_rcp_ps(_mm_rsqrt_ps(in));
        }
        FF_INLINE __m128 invSqrtFast(const __m128 in) {
            return _mm_rsqrt_ps(in);
        }
        FF_INLINE __m128 invPrecise(const __m128 in) {
            return _mm_div_ps(_mm_set1_ps(1.f), in);
        }
        FF_INLINE __m128 invSqrtPrecise(const __m128 in) {
            return _mm_div_ps(_mm_set1_ps(1.f), sqrtPrecise(in));
        }
        FF_INLINE __m128 sin(const __m128 in) {
            //return MathFun::sin_ps(in);
            return SseWrapper::sin(in).data;
        }
        FF_INLINE __m128 cos(const __m128 in) {
            //return MathFun::cos_ps(in);
            return SseWrapper::cos(in).data;
        }
        FF_INLINE __m128 log(const __m128 in) {
            return MathFun::log_ps(in);
        }
        FF_INLINE __m128 exp(const __m128 in) {
            return MathFun::exp_ps(in);
        }
        FF_INLINE void sincos(const __m128 in, __m128& outSin, __m128& outCos) {
            MathFun::sincos_ps(in, &outSin, &outCos);
        }    
        FF_INLINE __m128 acos(const __m128 in) {
            return SseWrapper::_mm_acos_ps(in);
        }
        FF_INLINE __m128 atan(const __m128 in) {
            return SseWrapper::atan(in).data;
        }    

        FF_INLINE __m128 powFast(const __m128 base, const float exponent) {
            return FastPow::fastPow(base, _mm_set1_ps(exponent));
        }
        FF_INLINE __m128 powFast(const __m128 base, const __m128 exponent) {
            return FastPow::fastPow(base, exponent);
        }

        FF_INLINE __m128 exp2fast(const __m128 x) {
            return FastPow::exp2f4(x); //TODO napsat tu skalarni verzi?
        }
        FF_INLINE __m128 log2fast(const __m128 x) {
            return FastPow::log2f4(x); //TODO napsat tu skalarni verzi?
        }



        FF_INLINE float sqrtFast(const float in) {
            __m128 temp = _mm_rcp_ss(_mm_rsqrt_ss(_mm_set1_ps(in)));
            return *(float*)&temp;
        }
        FF_INLINE float invFast(const float in) {
            __m128 temp = _mm_rcp_ss(_mm_set1_ps(in));
            return *(float*)&temp;
        }
        FF_INLINE float invSqrtFast(const float in) {
            __m128 temp = _mm_rsqrt_ss(_mm_set1_ps(in));
            return *(float*)&temp;
        }
        FF_INLINE float sqrtPrecise(const float in) {
            return ::sqrtf(in);
        }
        FF_INLINE float invSqrtPrecise(const float in) {
            return 1.f/sqrtf(in);
        }

        FF_INLINE void sincos(const float in, float& outSin, float& outCos) {
            __m128 temp1, temp2;
            MathFun::sincos_ps(_mm_set1_ps(in), &temp1, &temp2); //TODO napsat skalarni verzi tohohle
            outSin = *(float*)&temp1;
            outCos = *(float*)&temp2;   
        }
        FF_INLINE float atan2(const float y, const float x) {    //todo napsat skalar verzi?
            //__m128 temp = SseWrapper::atan2(_mm_set1_ps(x), _mm_set1_ps(y)).data;
            //return *(float*)&temp;
            return ::atan2(y, x);   //TODO namerit, pripadne vymenit
        }
        FF_INLINE __m128 atan2(const __m128 _y, const __m128 _x) {    //todo napsat skalar verzi?
            //__m128 temp = SseWrapper::atan2(_mm_set1_ps(x), _mm_set1_ps(y)).data;
            //return *(float*)&temp;
            const float* x = (float*)&_x;
            const float* y = (float*)&_y;
            return _mm_set_ps(
                ::atan2(y[3], x[3]),
                ::atan2(y[2], x[2]),
                ::atan2(y[1], x[1]),
                ::atan2(y[0], x[0])
                );
        }


        FF_INLINE float acos(const float x) {
            return ::acos(x);   //TODO namerit, pripadne vymenit
        }
        FF_INLINE float cos(const float x) {
            return ::cos(x);   //TODO namerit, pripadne vymenit
        }
        FF_INLINE float exp(const float x) {
            __m128 temp = exp(_mm_set1_ps(x));
            return *(float*)&temp;
        }

        FF_INLINE float importanceExp(const float _x) {
           //return fasterexp(_x);
            return fastexp(_x);
            //return ::expf(_x);
            //return (float&) MathFun::importanceExp(_mm_set1_ps(_x));
        }
        FF_INLINE __m128 importanceExp(const __m128 _x) {
            //return fasterexp(_x);
            return my_fastexp(_x);
            //return ::expf(_x);
            //return (float&) MathFun::importanceExp(_mm_set1_ps(_x));
        }

        FF_INLINE float sin(const float x) {
            return ::sin(x);   //TODO namerit, pripadne vymenit
        }
        FF_INLINE float log(const float x) {
            return ::log(x);   //TODO namerit, pripadne vymenit
        }
        FF_INLINE float tan(const float x) {
            return ::tan(x);   //TODO namerit, pripadne vymenit
        }
        FF_INLINE float atan(const float x) {
            return ::atan(x);   //TODO namerit, pripadne vymenit
        }
        FF_INLINE float powFast(const float base, const float exponent) {
            return *(const float*)&FastPow::fastPow(_mm_set1_ps(base), _mm_set1_ps(exponent));
            //TODO napsat tu skalarni verzi?
        }
        FF_INLINE float powPrecise(const float base, const float exponent) {
            return powf(base, exponent);
        }

        FF_INLINE float exp2fast(const float x) {
            return *(float*)&FastPow::exp2f4(_mm_set1_ps(x)); //TODO napsat tu skalarni verzi?
        }

        FF_INLINE float log2fast(const float x) {
            return *(float*)&FastPow::log2f4(_mm_set1_ps(x)); //TODO napsat tu skalarni verzi?
        }

        //TODo optimalizace - je to potreba nekde?

        /// \brief Calculates logarithm of float in given base
        /// \param number number from which the logarithm is computed
        /// \param base base of the logarithm
        /// \return log_base(number)
        FF_INLINE float log(const float number, const float base) {
            //TODO overit vzorec na par prikladech
            return ::log(number)/::log(base);  
        }

        template<int TExp>
        FF_INLINE float pow(const float value);

        template<>
        FF_INLINE float pow<2>(const float value) {
            return value*value;
        }

        template<>
        FF_INLINE float pow<3>(const float value) {
            return value*value*value;
        }

        template<>
        FF_INLINE float pow<4>(const float value) {
            return pow<2>(pow<2>(value));
        }
        template<>
        FF_INLINE float pow<5>(const float value) {
            return pow<2>(pow<2>(value))*value;
        }

        FF_INLINE float asinh(const float x) {
            return log(x + sqrt(x*x + 1));
        }
    }
}

#pragma warning(pop)

#endif