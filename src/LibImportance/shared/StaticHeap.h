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


#ifndef ___STATIC_HEAP__
#define ___STATIC_HEAP__

#include "StaticArray.h"

#pragma warning(push)
#pragma warning(disable:4127)

namespace Importance {

    /// \brief Default comparator for heap which uses element operator <
    template<class TItem>
    class OperatorLessComparator {
    public:
        /// \brief Returns true, if x < y, based on its overloaded operator<
        /// \param x first item to compare
        /// \param y second item to compare
        /// \return true iff x < y
        IMPORTANCE_INLINE bool operator()( const TItem& x, const TItem& y ) const {
            return x < y;
        }
    };


    /// \brief Priority queue implemented as heap in static (local) array, templated with maximal 
    ///        queue size, stored item and its comparator. The heap is min-heap when the 
    ///        comparator has the property that comp(x, y) == true iff x < y
    template <class TItem, class TStorage, class TComparator = OperatorLessComparator<TItem>>
    class StaticHeapBase {

        /// \brief Static storage of data for the heap
        TStorage data;

        /// \brief Comparator, used to compare heap elements. Comparator class must have 
        ///        overloaded operator() which takes 2 TItem instances as argument and returns 
        ///        true if first is strictly lower than second
        TComparator comp;   

        /// \brief Number of elements currently stored
        int stored;

    public:

        IMPORTANCE_INLINE StaticHeapBase(const int stored = 0) : stored(stored) { }

        IMPORTANCE_INLINE StaticHeapBase(const int stored, const TStorage& storage) : data(storage), stored(stored) { }

        /// \brief Sets new comparator for this heap
        /// \param comp new comparator for this heap
        IMPORTANCE_INLINE void setComparator(const TComparator& comp) {
            this->comp = comp;
        }

        /// \brief Returns number of elements currently stored in this heap
        /// \return Number of elements currently stored in this heap
        IMPORTANCE_INLINE int size() const {
            return stored;
        }

        /// \brief Inserts item into the queue
        /// \param item item to insert
        IMPORTANCE_INLINE void insert(const TItem& item) {
            IMPORTANCE_ASSERT(!isFull());
            data[stored] = item;
            int index = stored;
            ++stored;
            int tempIndex;
            while(index != 0 && comp(data[index], data[tempIndex = (index-1)>>1])) {
                std::swap(data[index], data[tempIndex]);
                index = tempIndex;
            }
        }

        /// \brief returns the "smallest" element from the queue
        /// \return Element from the top of the heap
        IMPORTANCE_INLINE TItem removeTop() {
            IMPORTANCE_ASSERT(!isEmpty());
            TItem retVal = data[0];
            data[0] = data[--stored];
            // repair heap

            int index = 0;
            while(true) {
                const int left = (index<<1) + 1;
                const int right = left+1;
                int largest = index;

                if(right < stored) {
                    if(comp(data[right], data[largest])) {
                        largest = right;
                    }
                    if(comp(data[left], data[largest])) {
                        largest = left;
                    }
                } else if(left < stored && comp(data[left], data[largest])) {
                    largest = left;
                }

                if(largest == index) {
                    break;
                }
                std::swap(data[largest], data[index]);
                index = largest;
            }
            return retVal;
        }

        /// \brief Replaces the element on top of this heap and repairs the heap, if needed
        /// \param item item to replace the top of heap with
        IMPORTANCE_INLINE void replaceTop(const TItem& item) {
            data[0] = item;
            // heap repairTop
            int index = 0;
            while(true) {
                const int left = (index<<1) + 1;
                const int right = left+1;
                int largest = index;


                if(right < stored) {
                    if(comp(data[right], data[largest])) {
                        largest = right;
                    }
                    if(comp(data[left], data[largest])) {
                        largest = left;
                    }
                } else if(left < stored && comp(data[left], data[largest])) {
                    largest = left;
                }


                if(largest == index) {
                    break;
                }
                std::swap(data[largest], data[index]);
                index = largest;        
            }
        }

        /// \brief Returns constant reference to the smallest item in the queue without 
        ///        dequeuing it
        /// \return constant reference to the smallest item in the queue 
        IMPORTANCE_INLINE const TItem& peek() const {
            IMPORTANCE_ASSERT(!isEmpty());
            return data[0];
        }

        /// \brief Returns reference to the smallest item in the queue without  dequeuing it
        /// \return reference to the smallest item in the queue 
        IMPORTANCE_INLINE TItem& peek() {
            IMPORTANCE_ASSERT(!isEmpty());
            return data[0];
        }

        /// \brief Tells if the queue is full (because it uses local array it cannot get resized)
        /// \return true if the queue is already full
        IMPORTANCE_INLINE bool isFull() const {
            return stored == data.size();
        }

        /// \brief Tells if the queue is empty
        /// \return true if size of this queue is 0
        IMPORTANCE_INLINE bool isEmpty() const {
            return stored == 0;
        }

        /// \brief Returns element at given position
        /// \param index index of element to return
        /// \return element at position index
        IMPORTANCE_INLINE const TItem& operator[](const unsigned int index) const {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(stored));
            return data[index];
        }
    };

    template<class TItem, int TSize, class TComparator = OperatorLessComparator<TItem>>
    class StaticHeap : public StaticHeapBase<TItem, ImStaticArray<TItem, TSize>, TComparator> {
    public:
        /// \brief Constructs an empty heap
        IMPORTANCE_INLINE StaticHeap() : StaticHeapBase(0) {
            static_assert(TSize >= 0, "Heap must have positive size");
        }

        /// \brief Constructs an empty heap with given comparator
        /// \param comp comparator to be used for this heap
        IMPORTANCE_INLINE StaticHeap(const TComparator& comp) : StaticHeapBase(0), comp(comp) {   
            static_assert(TSize >= 0, "Heap must have positive size");
        }
    };


    template<class TItem>
    class PtrAndSize{
    protected:
        TItem* ptr;
        int count;
    public:

        IMPORTANCE_INLINE PtrAndSize() { }

        IMPORTANCE_INLINE PtrAndSize(TItem* ptr, const int count) : ptr(ptr), count(count) {
            IMPORTANCE_ASSERT(count > 0);
        }

        IMPORTANCE_INLINE const TItem& operator[](const int index) const {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(count));
            return ptr[index];
        }
        IMPORTANCE_INLINE TItem& operator[](const int index) {
            IMPORTANCE_ASSERT(unsigned(index) < unsigned(count));
            return ptr[index];
        }

        IMPORTANCE_INLINE int size() const {
            return count;
        }
    };

    template<class TItem, class TComparator = OperatorLessComparator<TItem>>
    class ImExternalHeap : public StaticHeapBase<TItem, PtrAndSize<TItem>, TComparator> {
    public:
        IMPORTANCE_INLINE ImExternalHeap(TItem* data, const int heapSize, const int stored = 0) : StaticHeapBase(stored, PtrAndSize<TItem>(data, heapSize)) {   
            IMPORTANCE_ASSERT(heapSize >= stored);
        }
    };
}

#pragma warning(pop)
#endif