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
#include <vector>
#include "Stack.h"
/*<~API_TAG~>*/

#pragma warning(push)
#pragma warning(disable:4127)

namespace Importance {


    ///!!!!! column-wise

    /// \brief Class implementing 2D bitmap with templated pixel type
    template<typename TItem>
    class ImBitmap {
    protected:

        /// \brief internal data storage serialized to single array
        IStack<TItem> data;

        /// \brief width of the bitmap
        int width;

        /// \brief height of the bitmap
        int height;

        /// \brief 2D bitmap -> 1D array mapping function
        /// \param x x-coordinate to map
        /// \param y y-coordinate to map
        /// \return index into the serialized 1D array
        IMPORTANCE_INLINE int map(const int x, const int y) const {
            IMPORTANCE_ASSERT(unsigned(x) < unsigned(width) && unsigned(y) < unsigned(height));
            return y + x*height;
        }

    public:

        /// \brief Constructs empty bitmap (with zero width and height)
        ImBitmap() {
            this->width = 0;
            this->height = 0;
        }

        /// \brief Constructs bitmap with given height and width
        /// \param width width of the bitmap to construct
        /// \param height height of the bitmap to construct
        ImBitmap( const int width, const int height ) : data(width*height), width(width), height(height) {
            IMPORTANCE_ASSERT(width > 0); // height se zkontroluje v array
        }

        /// \brief Copy constructor, copies content from the other bitmap
        /// \param other bitmap to copy from
        ImBitmap( const ImBitmap& other ) : data(other.data), width(other.width), height(other.height) {}

        /// \brief Assignment operator, copies content of the other bitmap into this bitmap
        /// \param other bitmap to copy from
        /// \return self
        ImBitmap& operator=( const ImBitmap& other ) {
            if(&other == this) {
                IMPORTANCE_ASSERT(false);
                //return *this;
            }
            this->width = other.width;
            this->height = other.height;
            this->data = other.data;
            return *this;
        }

        void resize(const int width, const int height) {
            IMPORTANCE_ASSERT(width >= 0);    // height se zkontroluje v array
            if(this->width*this->height != width*height) {
                this->data.resize(width*height);
            }
            this->width = width;
            this->height = height;
        }

        /// \brief Returns non-const reference to given pixel, thus allowing modifications 
        ///        in bitmap
        /// \param x x-coordinate of pixel to return
        /// \param y y-coordinate of pixel to return
        /// \return non-const reference to pixel at coordinates [x, y]
        IMPORTANCE_INLINE TItem& get(const int x, const int y) {
            return data[map(x, y)];
        }

        /// \brief Returns const reference to given pixel
        /// \param x x-coordinate of pixel to return
        /// \param y y-coordinate of pixel to return
        /// \return const reference to pixel at coordinates [x, y]
        IMPORTANCE_INLINE const TItem& get(const int x, const int y) const {
            return data[map(x, y)];
        }

        /// \brief Returns non-const reference to given pixel, thus allowing modifications 
        ///        in bitmap. Equivalent to get
        /// \param x x-coordinate of pixel to return
        /// \param y y-coordinate of pixel to return
        /// \return non-const reference to pixel at coordinates [x, y]
        IMPORTANCE_INLINE TItem& operator()(const int x, const int y) {
            return get(x, y);
        }

        /// \brief Returns const reference to given pixel. Equivalent to get
        /// \param x x-coordinate of pixel to return
        /// \param y y-coordinate of pixel to return
        /// \return const reference to pixel at coordinates [x, y]
        IMPORTANCE_INLINE const TItem& operator()(const int x, const int y) const {
            return get(x, y);
        }
   
        /// \brief returns width of this bitmap
        /// \return width of this bitmap
        IMPORTANCE_INLINE int getWidth() const {
            return this->width;
        }
        
        /// \brief returns height of this bitmap
        /// \return height of this bitmap
        IMPORTANCE_INLINE int getHeight() const {
            return this->height;
        }

        IMPORTANCE_INLINE int numberOfElements() const {
            return width*height;
        }
        IMPORTANCE_INLINE TItem* ptr() {
            return &this->data[0];
        }
        IMPORTANCE_INLINE const TItem* ptr() const {
            return &this->data[0];
        }
        IMPORTANCE_INLINE TItem& getNth(const int index) {
            return this->data[index];
        }
        IMPORTANCE_INLINE const TItem& getNth(const int index) const {
            return this->data[index];
        }

        /// \brief returns total number of pixels in this bitmap
        /// \return number of pixels in this bitmap
        IMPORTANCE_INLINE int pixelCount() const {
            return this->width * this->height;
        }

        IMPORTANCE_INLINE void fill(const TItem& value) {
            for(int i = 0; i < data.size(); ++i) {
                data[i] = value;
            }
        }      
        IMPORTANCE_INLINE typename IStack<TItem>::const_iterator getDataIterator(const int x = 0, const int y = 0) const {
            return this->data.cbegin() + map(x, y);
        } 
    };
}

#pragma warning(pop)