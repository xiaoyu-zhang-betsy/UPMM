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


#include "../shared/Config.h"

#ifdef CUSTOM_RANDOM
/* 
   A C-program for MT19937-64 (2004/9/29 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/

#include "random.h"

namespace Importance {

    Random::Random() {
	    mti=MT_N+1; /* mti==N+1 means mt[N] is not initialized */
	    seed();
    }

    Random::Random(Random *random) {
	    mti=MT_N+1; /* mti==N+1 means mt[N] is not initialized */
	    seed(random);
    }

    Random::Random(uint64_t seedval) {
	    mti=MT_N+1; /* mti==N+1 means mt[N] is not initialized */
	    seed(seedval);
    }

    void Random::seed(uint64_t s) {
        mt[0] = s;
        for (mti=1; mti<MT_N; mti++) 
            mt[mti] =  (6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
    }

    void Random::seed(Random *random) {
	    uint64_t buf[MT_N];
	    for (int i=0; i<MT_N; ++i)
		    buf[i] = random->nextULong();
	    seed(buf, MT_N);
    }

    void Random::set(Random *random) {
	    for (int i=0; i<MT_N; ++i)
	        mt[i] = random->mt[i];
	    mti = random->mti;
    }

    void Random::seed(uint64_t *init_key, uint64_t key_length) {
        uint64_t i, j, k;
        seed(19650218ULL);
        i=1; j=0;
        k = (MT_N>key_length ? MT_N : key_length);
        for (; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 3935559000370003845ULL))
              + init_key[j] + j; /* non linear */
            i++; j++;
            if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
            if (j>=key_length) j=0;
        }
        for (k=MT_N-1; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * 2862933555777941757ULL))
              - i; /* non linear */
            i++;
            if (i>=MT_N) { mt[0] = mt[MT_N-1]; i=1; }
        }

        mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */ 
    }


    uint64_t Random::nextULong() {
        int i;
        uint64_t x;
        static uint64_t mag01[2]={0ULL, MT_MATRIX_A};

        if (mti >= MT_N) { /* generate MT_N words at one time */
            /* if init_genrand64() has not been called, */
            /* a default initial seed is used     */
            if (mti == MT_N+1) 
                seed(5489ULL); 

            for (i=0;i<MT_N-MT_M;i++) {
                x = (mt[i]&MT_UPPER_MASK)|(mt[i+1]&MT_LOWER_MASK);
                mt[i] = mt[i+MT_M] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
            }
            for (;i<MT_N-1;i++) {
                x = (mt[i]&MT_UPPER_MASK)|(mt[i+1]&MT_LOWER_MASK);
                mt[i] = mt[i+(MT_M-MT_N)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
            }
            x = (mt[MT_N-1]&MT_UPPER_MASK)|(mt[0]&MT_LOWER_MASK);
            mt[MT_N-1] = mt[MT_M-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];

            mti = 0;
        }
  
        x = mt[mti++];

        x ^= (x >> 29) & 0x5555555555555555ULL;
        x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
        x ^= (x << 37) & 0xFFF7EEE000000000ULL;
        x ^= (x >> 43);

        return x;
    }

    uint32_t Random::nextUInt(uint32_t n) {
	    /* Determine a bit mask */
	    uint32_t result, bitmask = n;
	    bitmask |= bitmask >> 1;
	    bitmask |= bitmask >> 2;
	    bitmask |= bitmask >> 4;
	    bitmask |= bitmask >> 8;
	    bitmask |= bitmask >> 16;

	    /* Generate numbers until one in [0, n) is found */
	    while ((result = (uint32_t) (nextULong() & bitmask)) >= n)
		    ;

	    return result;
    }

    size_t Random::nextSize(size_t n) {
	    /* Determine a bit mask */
	    size_t result, bitmask = n;
	    bitmask |= bitmask >> 1;
	    bitmask |= bitmask >> 2;
	    bitmask |= bitmask >> 4;
	    bitmask |= bitmask >> 8;
	    bitmask |= bitmask >> 16;

    #if defined(WIN64) || defined(__LINUX__) || defined(__OSX__)
	    if (sizeof(size_t) > 4)
		    bitmask |= bitmask >> 32;
    #endif

	    /* Generate numbers until one in [0, n) is found */
	    while ((result = (size_t) (nextULong() & bitmask)) >= n)
		    ;

	    return result;
    }

    #if defined(DOUBLE_PRECISION)
    Float Random::nextFloat() { 
	    return (Float) ((nextULong() >> 11) * (1.0/9007199254740992.0));
    }
    #else
    Float Random::nextFloat() {
	    /* Trick from MTGP: generate an uniformly distributed 
           single precision number in [1,2) and subtract 1. */
        union {
		    uint32_t u;
		    float f;
        } x;
        x.u = ((nextULong() & 0xFFFFFFFF) >> 9) | 0x3f800000UL;
        return x.f - 1.0f;
    }
    #endif

}

#endif