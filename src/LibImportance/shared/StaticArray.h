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


#ifndef ___STATIC_ARRAY___
#define ___STATIC_ARRAY___

#include <vector>
#include <xutility>
#include "Config.h"
#include "IIterator.h"

namespace Importance {

    template <typename Type, int TSize>
    class ImStaticArray  {
    protected: 

        Type data[TSize];

    public:

        //typedef IteratorBase<Type, StaticArray<Type, TSize>> IIterator;
        //typedef IteratorBase<const Type, const StaticArray<Type, TSize>> ConstIterator;

        IMPORTANCE_INLINE ImStaticArray() { }

        IMPORTANCE_INLINE explicit ImStaticArray(const Type& initValue) {
            for(int i = 0; i < TSize; ++i) {
                data[i] = initValue;
            }
        }

        /// \brief Deep copy copy constructor
        ImStaticArray(const ImStaticArray& second) {  
            for(int i = 0; i < TSize; ++i) {
                this->data[i] = second.data[i];
            }
        } 

        ImStaticArray& operator= (const ImStaticArray& second) {
            IMPORTANCE_ASSERT(this != &second);

            for(int i = 0; i < TSize; ++i) {
                this->data[i] = second.data[i];
            }
            return *this;
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
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(TSize));
            return data[index];
        }

        /// \brief Returns a element of the array at given index.
        /// \param index zero-based index of element to return
        /// \return reference to the element at index position
        IMPORTANCE_INLINE Type& get(const int index) {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(TSize));
            return data[index]; 
        }
    
        /// \brief returns size of this array
        /// \return size of this array
        IMPORTANCE_INLINE int size() const {
            return TSize;
        }

        /// \brief returns pointer to the beginning of internal continuous storage of elements
        /// \return pointer to the first element in which data are stored
        IMPORTANCE_INLINE Type* ptr() {
            return data;
        }

        /// \brief Returns constant pointer to the beginning of internal continuous storage 
        ///        of elements
        /// \return constant pointer to the first element in which data are stored
        IMPORTANCE_INLINE const Type* ptr() const {
            return data;
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

#endif