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
#include "Vector2.h"
#include "Vector3.h"
#include "../Shared/NumberUtils.h"
#include "Utils.h"

namespace Importance {
    template <int M, int N, typename T> struct Matrix {
    public:
        T m[M][N];

        /// Indexing operator
        inline T &operator()(int i, int j) { return m[i][j]; }

        /// Indexing operator (const verions)
        inline const T & operator()(int i, int j) const { return m[i][j]; }    
    };

    class Matrix3x3 : public Matrix<3,3,Float> {
    public:
        /// Matrix-vector multiplication
        inline Vector3 operator*(const Vector3 &v) const {
            return Vector3(
                m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
        }
    };

	// Symetric, positive-definite matrix 2x2
	template< typename T >
	struct SPDMatrix2x2
	{
		T m[ 3 ];

		IMPORTANCE_INLINE SPDMatrix2x2()
		{
			m[ 0 ] = m[ 1 ] = m[ 2 ] = T(0.0f);
		}

		IMPORTANCE_INLINE SPDMatrix2x2( const T & m_11, const T & m_12, const T & m_22 )
		{
			m[ 0 ] = m_11;
			m[ 1 ] = m_12;
			m[ 2 ] = m_22;
		}
		
		// Sufficient statistics
		template < typename T2 >
		IMPORTANCE_INLINE SPDMatrix2x2( const Vector2t< T2 > & vector )
		{
			m[ 0 ] = vector.x * vector.x;
			m[ 1 ] = vector.x * vector.y;
			m[ 2 ] = vector.y * vector.y;
		}

		// Indexing operator
		IMPORTANCE_INLINE T &operator()(int i, int j) 
		{
			return m[ i + j ];
		}

		// Indexing operator (const)
		IMPORTANCE_INLINE const T &operator()(int i, int j) const
		{
			return m[ i + j ];
		}

		// Matrix += Matrix
		IMPORTANCE_INLINE SPDMatrix2x2 operator+=( const SPDMatrix2x2 & matrix )
		{
			m[ 0 ] += matrix.m[ 0 ];
			m[ 1 ] += matrix.m[ 1 ];
			m[ 2 ] += matrix.m[ 2 ];
			return *this;
		}

		// Matrix *= Scalar
		IMPORTANCE_INLINE SPDMatrix2x2 operator*=( const T & scalar )
		{
			m[ 0 ] *= scalar;
			m[ 1 ] *= scalar;
			m[ 2 ] *= scalar;
			return *this;
		}

		// Matrix /= Scalar
		IMPORTANCE_INLINE SPDMatrix2x2 operator/=( const T & scalar )
		{
			m[ 0 ] /= scalar;
			m[ 1 ] /= scalar;
			m[ 2 ] /= scalar;
			return *this;
		}

		// Matrix + Matrix
		IMPORTANCE_INLINE SPDMatrix2x2 operator+( const SPDMatrix2x2 & matrix ) const
		{
			SPDMatrix2x2 mat;
			mat.m[ 0 ] = m[ 0 ] + matrix.m[ 0 ];
			mat.m[ 1 ] = m[ 1 ] + matrix.m[ 1 ];
			mat.m[ 2 ] = m[ 2 ] + matrix.m[ 2 ];
			return mat;
		}



		// Matrix * Scalar
		IMPORTANCE_INLINE SPDMatrix2x2 operator*( const T & scalar ) const
		{
			SPDMatrix2x2 mat;
			mat.m[ 0 ] = m[ 0 ] * scalar;
			mat.m[ 1 ] = m[ 1 ] * scalar;
			mat.m[ 2 ] = m[ 2 ] * scalar;
			return mat;
		}

		// Matrix / Scalar
		IMPORTANCE_INLINE SPDMatrix2x2 operator/( const T & scalar ) const
		{
			SPDMatrix2x2 mat;
			mat.m[ 0 ] = m[ 0 ] / scalar;
			mat.m[ 1 ] = m[ 1 ] / scalar;
			mat.m[ 2 ] = m[ 2 ] / scalar;
			return mat;
		}

		// Matrix * Vector
		template < typename T2 >
		IMPORTANCE_INLINE Vector2t< T2 > operator*( const Vector2t< T2 > & vector ) const
		{
			return Vector2t< T2 >( m[ 0 ] * vector.x + m[ 1 ] * vector.y, m[ 1 ] * vector.x + m[ 2 ] * vector.y );
		}

		// Quadratic form: vector^T * Matrix * vector
		template < typename T2 >
		IMPORTANCE_INLINE T quadraticForm( const Vector2t< T2 > & vector ) const
		{
			const T tmp = vector.x * vector.y * m[ 1 ];
			return vector.x * vector.x * m[ 0 ] + tmp + tmp + vector.y * vector.y * m[ 2 ];
		}

		// Determinant
		IMPORTANCE_INLINE T determinant() const
		{
			return m[ 0 ] * m[ 2 ] - m[ 1 ] * m[ 1 ];
		}

        // Trace
        IMPORTANCE_INLINE T trace() const {
            return m[ 0 ] + m[ 2 ];
        }

		// Invert matrix
		IMPORTANCE_INLINE void invert()
		{
			const T inv_det = T(1.0f) / determinant();
			const T tmp = m[ 0 ];
			m[ 0 ] = m[ 2 ] * inv_det;
			m[ 1 ] *= -inv_det;
			m[ 2 ] = tmp * inv_det;
		}
		
        IMPORTANCE_INLINE SPDMatrix2x2 getInverse() const {
            const T inv_det = T(1.0f) / determinant();
            return SPDMatrix2x2(m[ 2 ] * inv_det, m[ 1 ] * -inv_det, m[ 0 ] * inv_det);
        }


		// Turn matrix to scaled identity matrix
		IMPORTANCE_INLINE void identity( const T & scale = T( 1 ) )
		{
			m[ 0 ] = m[ 2 ] = T( scale );
			m[ 1 ] = T( 0.0f );
		}

        IMPORTANCE_INLINE bool isReal() const {
            return Importance::isReal(m[0]) && Importance::isReal(m[1]) && Importance::isReal(m[2]);
        }

        static IMPORTANCE_INLINE SPDMatrix2x2 lerp(const SPDMatrix2x2& a, const SPDMatrix2x2& b, const float amountB) {
            return SPDMatrix2x2(
                Importance::lerp(a.m[0], b.m[0], amountB),
                Importance::lerp(a.m[1], b.m[1], amountB),
                Importance::lerp(a.m[2], b.m[2], amountB)
                );
        }

	};

    template< typename T >
    IMPORTANCE_INLINE SPDMatrix2x2<T> operator*( const T & scalar, const SPDMatrix2x2<T> & m )
    {
        return m * scalar;
    }

}