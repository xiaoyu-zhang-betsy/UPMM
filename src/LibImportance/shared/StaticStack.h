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


#ifndef ___STATIC_STACK___
#define ___STATIC_STACK___

#include "StaticArray.h"

namespace Importance {

    /// \brief Minimalistic, non-resizable templated stack without heap memory allocation. Has 
    ///        safety checks only in debug mode (via macros and assert())
    template <typename Type, int Tcapacity>
    class StaticStack {
    protected: 

        /// \brief internal storage of data
        ImStaticArray<Type, Tcapacity> storage;
    
        /// \brief Number of currently stored elements
        int _size;

    public:

        /// \brief constructs an empty stack
        IMPORTANCE_INLINE StaticStack() :_size(0) {
            static_assert(Tcapacity >= 0, "Stack must have positive size");
        }

        /// \brief pushes one element to the top of the stack
        /// \param item item to store
        IMPORTANCE_INLINE void push(const Type item) {
            IMPORTANCE_ASSERT(!isFull());
            storage[_size] = item;
            ++_size;
        }

        /// \brief Removes and returns one element from stack
        /// \return Element from the top of the stack
        IMPORTANCE_INLINE Type pop(){
            IMPORTANCE_ASSERT(!isEmpty());
            --_size;
            return storage[_size];
        }

        IMPORTANCE_INLINE const Type* ptr() const {
            return storage.ptr();
        }

        /// \brief Returns constant reference to element from the top of the stack without 
        ///        removing it
        /// \return constant reference to top element of the stack
        IMPORTANCE_INLINE const Type& peek() const {
            IMPORTANCE_ASSERT(!isEmpty());
            return storage[_size-1];
        }

        /// \brief Returns non-const reference to element from the top of the stack without 
        ///        removing it
        /// \return reference to top element of the stack
        IMPORTANCE_INLINE Type& peek() {
            IMPORTANCE_ASSERT(!isEmpty());
            return storage[_size-1];
        }

        /// \brief tells if the stack is empty
        /// \return true if stack holds 0 elements, false otherwise
        IMPORTANCE_INLINE bool isEmpty() const{
            return _size==0;
        }

        IMPORTANCE_INLINE bool isFull() const{
            return _size==Tcapacity;
        }

        IMPORTANCE_INLINE int size() const {
            return this->_size;
        }

        IMPORTANCE_INLINE void resize(const int newSize) {
            IMPORTANCE_ASSERT(unsigned(newSize) <= Tcapacity);
            this->_size = newSize;
        }

        /// \brief Returns a constant element of the array at given index. Equivalent to the 
        ///        method get
        /// \param index zero-based index of element to return
        /// \return constant reference to the element at index position
        IMPORTANCE_INLINE const Type& operator[](const int index) const {
            return this->get(index);
        }

        /// \brief Returns a constant element of the array at given index. Equivalent to the 
        ///        method get
        /// \param index zero-based index of element to return
        /// \return reference to the element at index position
        IMPORTANCE_INLINE Type& operator[](const int index) {
            return this->get(index);
        }

        /// \brief Returns a constant element of the array at given index.
        /// \param index zero-based index of element to return
        /// \return constant reference to the element at index position
        IMPORTANCE_INLINE const Type& get(const int index) const {
            IMPORTANCE_ASSERT(index < this->_size);
            return storage[index];
        }

        /// \brief Returns a element of the array at given index.
        /// \param index zero-based index of element to return
        /// \return reference to the element at index position
        IMPORTANCE_INLINE Type& get(const int index) {
            IMPORTANCE_ASSERT(index < this->_size);
            return storage[index];
        }
    };
}

#endif