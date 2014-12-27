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
#include "Config.h"
/*<~API_TAG~>*/


namespace Importance {

/// \brief IIterator for the array, allowing to use it with STL algorithms
template<class IterType>
class IIterator {
protected:

    IterType* ptr;

#ifdef LIBIMP_DEBUG    
    IterType* begin;
    IterType* end;

    IMPORTANCE_INLINE void assertSameArray(const IIterator& other) const {
        IMPORTANCE_ASSERT(begin == other.begin);
        IMPORTANCE_ASSERT(end == other.end);
    }
#else
    IMPORTANCE_INLINE void assertSameArray(const IIterator&) const { }
#endif

public:

    template<class TArray>
    IMPORTANCE_INLINE IIterator(TArray* array, const size_t ptr) {
        this->ptr = array->ptr()+ ptr;
#ifdef LIBIMP_DEBUG
        this->begin = array->ptr();
        this->end = array->ptr() + array->size();
#endif
    }
    IMPORTANCE_INLINE explicit IIterator(IterType* ptr) {
        this->ptr = ptr;
#ifdef LIBIMP_DEBUG
        this->begin = ptr;
        this->end = ptr+1;
#endif
    }

    IMPORTANCE_INLINE IIterator(IterType* ptr, const size_t offset, const size_t arraySize) {
        this->ptr = ptr+offset;
#ifdef LIBIMP_DEBUG
        this->begin = ptr;
        this->end = ptr+arraySize;
#endif
    }

    IMPORTANCE_INLINE IIterator() { }

    IMPORTANCE_INLINE operator IIterator<const IterType>() const {
        return *(IIterator<const IterType>*)this;
    }

    IMPORTANCE_INLINE bool operator==(const IIterator& other) const {
        assertSameArray(other);
        return ptr == other.ptr;
    }
    IMPORTANCE_INLINE bool operator!=(const IIterator& other) const {
        return !((*this) == other);
    }
    IMPORTANCE_INLINE bool operator<(const IIterator& second) const {
        assertSameArray(second);
        return ptr < second.ptr;
    }
    IMPORTANCE_INLINE bool operator<=(const IIterator& second) const {
        assertSameArray(second);
        return ptr <= second.ptr;
    }
    IMPORTANCE_INLINE size_t operator-(const IIterator& second) const {
        assertSameArray(second);
        return ptr - second.ptr;
    }
    IMPORTANCE_INLINE IterType& operator*() {
        IMPORTANCE_ASSERT(ptr >= begin && ptr < end);
        return *ptr;
    }
    IMPORTANCE_INLINE const IterType& operator*() const {
        IMPORTANCE_ASSERT(ptr >= begin && ptr < end);
        return *ptr;
    }
    IMPORTANCE_INLINE IterType& operator[](int offset) {
        IMPORTANCE_ASSERT(ptr+offset >= begin && ptr+offset < end);
        return ptr[offset];
    }
    IMPORTANCE_INLINE const IterType& operator[](int offset) const {
        IMPORTANCE_ASSERT(ptr+offset >= begin && ptr+offset < end);
        return ptr[offset];
    }
    IMPORTANCE_INLINE IterType* operator->() {	
        IMPORTANCE_ASSERT(ptr >= begin && ptr < end);
        return (&**this);
    }
    IMPORTANCE_INLINE IIterator& operator++() {
        ++ptr;
        return *this;
    }
    IMPORTANCE_INLINE IIterator operator++(int) { 
        IIterator temp = *this;
        ++ptr;
        return temp;
    }
    IMPORTANCE_INLINE IIterator& operator--() {
        --ptr;
        return *this;
    }
    IMPORTANCE_INLINE IIterator operator--(int) { 
        IIterator temp = *this;
        --ptr;
        return temp;
    }
    IMPORTANCE_INLINE IIterator& operator+=(const size_t x) {
        ptr += x;
        return *this;
    }
    IMPORTANCE_INLINE IIterator operator+(const size_t x) const {
        IIterator res(*this);
        res.ptr += x;
        return res;
    }
    IMPORTANCE_INLINE IIterator operator-(const size_t x) const {
        IIterator res(*this);
        res.ptr -= x;
        return res;
    }

    // stl black magic
    typedef std::random_access_iterator_tag iterator_category;
    typedef IterType value_type;
    typedef size_t difference_type;
    typedef IterType* pointer;
    typedef IterType& reference;
};

}