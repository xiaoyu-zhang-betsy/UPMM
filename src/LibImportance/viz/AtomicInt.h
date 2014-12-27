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
#include <intrin.h>

namespace Importance {

    typedef __int64 int64;

//vraci result of operation
inline int64 atomicDecrement(volatile int64* base) {
    return _InterlockedDecrement64(base);
}

//vraci result of operation
inline int64 atomicIncrement(volatile int64* base) {
    return _InterlockedIncrement64(base);
}

//vraci result of operation
inline int64 atomicAdd(volatile int64* base, const int64 amount) {
    return _InterlockedExchangeAdd64(base, amount) + amount;
}

// opsano z embree
template<class T>
class AtomicInt {
protected: 
    volatile T data;

public:
    inline AtomicInt() {}
    inline AtomicInt(const T data) : data(data) { }
    inline AtomicInt& operator=(const T data) {
        this->data = data;
        return *this;
    }
    inline operator T() const { 
        return this->data; 
    }

public:

    inline T operator++() {
        return atomicIncrement(&data);
    }
    inline T operator--() {
        return atomicDecrement(&data);
    }
    inline T operator++(int) {
        return atomicIncrement(&data)-1;
    }
    inline T operator--(int) {
        return atomicDecrement(&data)+1;
    }
    inline T operator+=(const T amount) {
        return atomicAdd(&data, amount);
    }
};

//typedef AtomicInt<int> AtomicInt32;
typedef AtomicInt<int64> AtomicInt64;

}