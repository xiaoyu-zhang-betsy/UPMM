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

#include "Array.h"

namespace Importance {

    /// \brief Simple, fast templated stack, implemented using Corona::Array
    template <class Type, typename TInt = size_t>
    class IStack  {
        template<class Type> friend class IIterator;
    protected: 

        /// \brief internal storage of data
        IArray<Type> data;

        /// \brief current number of stored elements
        TInt stored;

        IMPORTANCE_INLINE Type* ptr() {
            return this->data.ptr();
        }
        IMPORTANCE_INLINE const Type* ptr() const {
            return this->data.ptr();
        }

    public:

        typedef IIterator<Type> iterator;
        typedef IIterator<const Type> const_iterator;
        
        /// \brief Creates new instance with empty storage
        IMPORTANCE_INLINE IStack() : stored(0) { }

        /// \brief A copy constructor, creates new instance and copies content of given stack 
        ///        into it
        /// \param second stack to copy
        IStack(const IStack& second) : data(second.data), stored(second.stored) { }

        /// \brief Created a stack with 0 elements, but preallocates space in the internal 
        ///        storage
        /// \param preallocate How many elements to preallocate
        explicit IStack(const size_t size) : stored(size), data(size) { }

        explicit IStack(const size_t size, const Type& val) : stored(size), data(size) {
            for(size_t i = 0; i < size; ++i) {
                data[i] = val;
            }
        }



        /// \brief An assignment operator, performs deep copy of the stack
        /// \param second stack from which to copy into this one
        /// \return self
        IStack& operator= (const IStack& second) {
            IMPORTANCE_ASSERT(this != &second);

            this->stored = second.stored;
            this->data = second.data;
            return *this;
        }

        /// \brief Returns a constant element of the stack at given index. Equivalent 
        ///        to the method get
        /// \param index zero-based index of element to return
        /// \return constant reference to the element at index position
        IMPORTANCE_INLINE const Type& operator[](const size_t index) const {
            return this->get(index);
        }

        /// \brief Returns an element of the stack at given index. Equivalent to the method get
        /// \param index zero-based index of element to return
        /// \return reference to the element at index position
        IMPORTANCE_INLINE Type& operator[](const size_t index) {
            return this->get(index);
        }

        /// \brief Returns a constant element of the stack at given index. 
        /// \param index zero-based index of element to return
        /// \return constant reference to the element at index position
        IMPORTANCE_INLINE const Type& get(const size_t index) const {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(this->stored));
            return this->data[index];
        }

        /// \brief Returns an element of the stack at given index. 
        /// \param index zero-based index of element to return
        /// \return reference to the element at index position
        IMPORTANCE_INLINE Type& get(const size_t index) {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(this->stored));
            return this->data[index];
        }

        /// \brief clears all entries (allocated storage is preserved)
        IMPORTANCE_INLINE void clear() {
            stored = 0;
        }

        /// \brief clears all entries as well as deletes allocated memory
        void dealloc() {
            stored = 0;
            data.dealloc();
        }

        /// \brief Tells if this stack is empty
        /// \return true if there are no elements stored in this stack
        IMPORTANCE_INLINE bool isEmpty() const {
            return size() == 0;
        }

        /// \brief Returns number of elements stored
        /// \return number of elements in the stack
        IMPORTANCE_INLINE TInt size() const {
            return stored;
        }

        /// \brief Pushes single element onto the stack. Reallocates internal storage if needed 
        ///        (with O(1) amortized time complexity)
        /// \param item item to push onto the stack
        IMPORTANCE_INLINE void push(const Type& item = Type()) {
            this->reserve(stored+1);
            this->data[int(stored++)] = item;
        }

        /// \brief returns reference to the top element in the stack
        /// \return reference to element at position (stored-1)
        IMPORTANCE_INLINE Type& peek() {
            return this->data[int(stored-1)];
        }

        /// \brief returns constant reference to the top element in the stack
        /// \return constant reference to element at position (stored-1)
        IMPORTANCE_INLINE const Type& peek() const {
            return this->data[int(stored-1)];
        }

        IMPORTANCE_INLINE void pop() {
            IMPORTANCE_ASSERT(!this->isEmpty());
            --stored;
        }

        /// \brief reserves at least given minimal capacity of the stack
        /// \param amount new minimal capacity of the stack
        IMPORTANCE_INLINE void reserve(const TInt amount) {
            IMPORTANCE_ASSERT(amount >= 0);
            if(this->data.size() < amount) {
                this->data.resize(int(std::max(amount, TInt(2*data.size()))));
            }
        }

        IMPORTANCE_INLINE void remove(const int index) {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(this->size()));
            for(int i = index+1; i < stored; ++i) {
                this->data[i-1] = this->data[i];
            }
            --this->stored;
        }   

        /// \brief reallocates the stack so that there is no unused space (number of elements 
        ///        stored is same as internal array capacity). Note that during execution 
        ///        of this method memory usage can actually spike because of the reallocation
        IMPORTANCE_INLINE void trim() {
            this->data.resize(stored);
        }

        /// \brief Pushes all elements from second stack onto this stack. Content of current 
        ///        stack remains unchanged at the bottom, content of second stack is copied 
        ///        after that so that so second stack top becomes new top
        /// \param other stack to copy onto this one
        IMPORTANCE_INLINE void pushAll(const IStack& other) {
            for(int i = 0; i < other.size(); ++i) {
                this->push(other[i]);
            }
        }

        IMPORTANCE_INLINE void pushAll(const IArray<Type>& other) {
            for(int i = 0; i < other.size(); ++i) {
                this->push(other[i]);
            }
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

        /// \brief Saves elements stored in this stack into an array. Overwrites previous content
        ///        of the array
        /// \param array Where to save
        void saveToArray(IArray<Type>& array) const {
            array.resize(this->stored);
            for(int i =0; i < this->stored; ++i) {
                array[i] = this->data[i];
            }
        }

        /// \brief sets number of elements stored. When newSize is bigger than current number 
        ///        of stored elements, value of new elements is undefined
        /// \param newSize new number of elements stored in the stack
        IMPORTANCE_INLINE void resize(const size_t newSize, const Type& initVal) {
            IMPORTANCE_ASSERT(newSize >= 0);
            const size_t oldSize = this->size();
            this->reserve(newSize);
            this->stored = newSize;
            for(size_t i = oldSize; i < newSize; ++i) {
                this->data[i] = initVal;
            }
        }
        IMPORTANCE_INLINE void resize(const size_t newSize) {
            IMPORTANCE_ASSERT(newSize >= 0);
            this->reserve(newSize);
            this->stored = newSize;
        }

        /// \brief removes an element from the middle of the stack in O(n) time
        /// \param index index of the element to remove
        IMPORTANCE_INLINE void remove(const size_t index) {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(this->size()));
            for(int i = index+1; i < stored; ++i) {
                this->data[i-1] = this->data[i];
            }
            --this->stored;
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