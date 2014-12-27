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

#include "../shared/Config.h"
#include <sstream>

namespace Importance {

    template<typename T>
    class Vector2t {
    public:
        typedef typename T::Scalar Scalar;
        typedef Vector2t<T> PointType;
        typedef Vector2t<T> VectorType;
        /// Number of dimensions
        const static int dim = 2;

        Scalar x, y;

        IMPORTANCE_INLINE Vector2t() {
#ifdef LIBIMP_DEBUG
        x = y = T::initVal;
#endif
        }

        IMPORTANCE_INLINE Vector2t(const Scalar x, const Scalar y) : x(x), y(y) { }

        IMPORTANCE_INLINE explicit Vector2t( Scalar val ) : x(val), y(val) {}

        /// Index into the vector's components
        Scalar &operator[](int i) {
            return (&x)[i];
        }

        /// Index into the vector's components (const version)
        Scalar operator[](int i) const {
            return (&x)[i];
        }

        Vector2t operator-(const Vector2t & a) const {
            return Vector2t( x - a.x, y - a.y );
        }
        Vector2t operator*(const Vector2t & a) const {
            return Vector2t( x * a.x, y * a.y );
        }

        Vector2t operator/(const Vector2t & a) const {
            return Vector2t( x / a.x, y / a.y );
        }

        Vector2t operator*( const Scalar & f ) const {
            return Vector2t( f * x, f * y );
        }

        Vector2t operator/( const Scalar & f ) const {
            const Scalar inverse = Scalar(1.f)/f;
            return Vector2t( x * inverse, y * inverse );
        }

        Vector2t operator+( const Vector2t & a ) const {
            return Vector2t( x + a.x, y + a.y );
        }

        Vector2t & operator*=( const Scalar & f ) {
            x *= f; y *= f;
            return *this;
        }

        Vector2t & operator/=( const Scalar & f ) {
            Scalar recip = (Scalar) 1 / f;
            x *= recip; y *= recip;
            return *this;
        }

        std::string toString() const {
            std::ostringstream str;
            str << "Vec[" << x << "," << y << "]";
            return str.str();
        }

        Scalar size() const {
            return std::sqrt(square());
        }
        Scalar square() const {
            return x*x + y*y;
        }

        IMPORTANCE_INLINE bool isValid() const {
            return isReal(x) && isReal(y);
        }
    };

    template<class T>
    IMPORTANCE_INLINE float dot(const Vector2t<T> a, const Vector2t<T> b) {
        return a.x*b.x + a.y*b.y;
    }

    template<typename T>
    IMPORTANCE_INLINE Vector2t<T> operator*( typename T::Scalar scalar, const Vector2t<T> & vector ) {
        return vector * scalar;
    }

    class VectorFloatTraits {
    public:
        typedef Float Scalar;
        static const Float initVal;
    };

    class VectorIntTraits {
    public:
        typedef int Scalar;
        static const int initVal = 0;
    };
    
	template< typename T >
	class VectorTraits {
		public:
			typedef T Scalar;
			static const T initVal;
	};

    typedef Vector2t<VectorFloatTraits> Vector2;
    typedef Vector2t<VectorIntTraits> Vector2i;
}