#pragma once

__m128 log2f4(const __m128 x);

__m128 exp2f4(__m128 x);

__forceinline __m128 fastPow(const __m128 x, const __m128 exponent) {
    return exp2f4(_mm_mul_ps(exponent, log2f4(x)));
}