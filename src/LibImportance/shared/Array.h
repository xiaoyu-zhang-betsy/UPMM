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

#include <xutility>
#include "Config.h"
#include "IIterator.h"

namespace Importance {

template <typename Type>
class IArray {
    template<class Type> friend class IIterator;
    template<class Type, class TInt> friend class IStack;
protected: 

    /// \brief Allocated size of data array
    size_t capacity;

    /// \brief Pointer to data stored in C array
    Type* data;

    IMPORTANCE_INLINE Type* ptr() {
        return data;
    }
    IMPORTANCE_INLINE const Type* ptr() const {
        return data;
    }

public:


    /// \brief Default constructor, creates array with zero size. No memory is allocated.
    IMPORTANCE_INLINE IArray() : capacity(0), data(NULL) { }

    /// \brief Deep copy copy constructor
    IArray(const IArray& second) {
        this->capacity = second.capacity;
        if(this->capacity > 0) {
            this->data = new Type[capacity];
            for(int i = 0; i < capacity; ++i) {
                this->data[i] = second.data[i];
            }
        } else {
            this->data = NULL;
        }
        
    }

    /// \brief Assignment operator, performs deep copy of second array to this array.
    /// \param second Array which to duplicate  
    /// \return reference to this array
    IArray& operator= (const IArray& second) {
        IMPORTANCE_ASSERT(capacity == 0 || data);
        IMPORTANCE_ASSERT(this != &second);
        if(this->capacity != second.capacity) {
            if(this->data) {
                delete [] this->data;
            }
            this->capacity = second.capacity;
            if(this->capacity == 0) {
                this->data = NULL;
            } else {
                this->data = new Type[capacity];
            }
            IMPORTANCE_ASSERT(capacity == 0 || data);
        }
        IMPORTANCE_ASSERT(capacity == 0 || data);
        for(int i = 0; i < capacity; ++i) {
            this->data[i] = second.data[i];
        }
        return *this;
    }

    /// \brief Constructs new array with given capacity
    /// \param capacity size of the array to be constructed
    explicit IArray(const size_t capacity) {
        IMPORTANCE_ASSERT(capacity >= 0);
        this->capacity = capacity;
        if(this->capacity > 0) {
            this->data = new Type[capacity];
        } else {
            this->data = NULL;
        }
    }


    ~IArray() {
        if(this->data) {
            delete[] this->data;
        }
    }

    /// \brief Returns a constant element of the array at given index. Equivalent to the 
    ///        method get
    /// \param index zero-based index of element to return
    /// \return constant reference to the element at index position
    IMPORTANCE_INLINE const Type& operator[](const size_t index) const {
        return this->get(index);
    }

    /// \brief Returns a constant element of the array at given index. Equivalent to the 
    ///        method get
    /// \param index zero-based index of element to return
    /// \return reference to the element at index position
    IMPORTANCE_INLINE Type& operator[](const size_t index) {
        return this->get(index);
    }

    /// \brief Returns a constant element of the array at given index.
    /// \param index zero-based index of element to return
    /// \return constant reference to the element at index position
    IMPORTANCE_INLINE const Type& get(const size_t index) const {
        IMPORTANCE_ASSERT(unsigned(index) < unsigned(size()));
        return data[index];
    }

    /// \brief Returns a element of the array at given index.
    /// \param index zero-based index of element to return
    /// \return reference to the element at index position
    IMPORTANCE_INLINE Type& get(const size_t index) {
        IMPORTANCE_ASSERT(unsigned(index) < unsigned(this->size()));
        return data[index];
    }

    /// \brief Resizes the array to exact new size. Elements 0 - min(new size, old size) 
    ///        remain unchanged.
    /// \param newSize new number of elements in the array
    void resize(const size_t newSize) {
		IMPORTANCE_ASSERT(newSize >= 0);
		if(newSize == this->capacity) {
			return;
		}

        if(newSize == 0) {
            dealloc();
        } else {
            Type* newData = new Type[newSize];
            size_t limit = std::min(newSize, capacity);

            for(size_t i = 0; i < limit; ++i) {
                newData[i] = this->data[i];
            }

            if(this->data != NULL) {
                delete [] this->data;   
            }
            this->data = newData;
            this->capacity = newSize;
        }
    }

    /// \brief deallocates internal storage, freeing memory, and setting size to 0.
    void dealloc() {
        if(this->data) {
            delete [] this->data;
            this->data = NULL;
            this->capacity = 0;
        } else {
            IMPORTANCE_ASSERT(capacity == 0);
        }
        
    }

    /// \brief returns size of this array
    /// \return size of this array
    IMPORTANCE_INLINE size_t size() const {
        return capacity;
    }

    /// \brief Returns constant iterator pointing at the first element of the array
    IMPORTANCE_INLINE IIterator<const Type> cbegin() const {
        return IIterator<const Type>(this, 0);
    }

    /// \brief Returns constant iterator pointing past the last element of the array
    IMPORTANCE_INLINE IIterator<const Type> cend() const {
        return IIterator<const Type>(this, size());
    }

    /// \brief Returns iterator pointing at the first element of the array
    IMPORTANCE_INLINE IIterator<Type> begin() {
        return IIterator<Type>(this, 0);
    }

    /// \brief Returns iterator pointing past the last element of the array
    IMPORTANCE_INLINE IIterator<Type> end() {
        return IIterator<Type>(this, size());
    }


    IMPORTANCE_INLINE const Type& back() const {
        return (*this)[size()-1];
    }
    IMPORTANCE_INLINE const Type& front() const {
        return (*this)[0];
    }
    IMPORTANCE_INLINE Type& back() {
        return (*this)[size()-1];
    }
    IMPORTANCE_INLINE Type& front() {
        return (*this)[0];
    }
};

}