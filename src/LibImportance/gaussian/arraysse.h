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


#ifndef __ARRAY_GAUSSIAN_SSE_H__
#define __ARRAY_GAUSSIAN_SSE_H__

#include <vector>
#include "../shared/Particle.h"
#include "../shared/frame.h"
#include "../shared/Utils.h"
#include "../shared/QuadVector.h"
#include "../shared/Vector2.h"
#include <malloc.h>

namespace Importance {

///-------------------------------------------------------------------------------------------------
/// Array of samples  
///-------------------------------------------------------------------------------------------------
template<class TScalar, int DIMENSION, typename TDistribution>
class ArrayGaussianSSE
{
    void operator=(const ArrayGaussianSSE&) { }
public:
	typedef Vector2t<VectorTraits< TScalar > > QuadVector2;

	IMPORTANCE_INLINE explicit ArrayGaussianSSE( const Frame & localFrame,  const Importance::Hit& hit, size_t reserve ):
		m_localFrame( localFrame ), m_hit( hit ), weight( 0 ), denom( 0 ), m_furthest( 0.0f )
	{
		m_samples.reserve( reserve );
	}

	IMPORTANCE_INLINE ~ArrayGaussianSSE()
	{
		_aligned_free( weight );
		_aligned_free( denom );
	}

	struct Sample
	{
		QuadVector2 pos, shifted_pos;
		TScalar weight;
		TScalar * gamma;
        TScalar * gamma_w;
		Vector3 original; // For debug purpose
		Vector3 photonStart;
        Float distance;

        //////////////////////////////////////////////////////////////////////////
        /// Adaption for stepwise E-M interface
        //////////////////////////////////////////////////////////////////////////

        IMPORTANCE_INLINE const QuadVector2 & getPoint() const { return pos; }

        IMPORTANCE_INLINE const TScalar & getValue() const { return weight; } 

        IMPORTANCE_INLINE Float getSingleValue() const { return weight[ 0 ]; }

        IMPORTANCE_INLINE void setValue( Float val ){ weight = TScalar( val ); }
	};  

    class WeightComparator {        
    public:
        IMPORTANCE_INLINE bool operator() (const Sample& l, const Sample& r) const {
            return l.weight[ 0 ] < r.weight[ 0 ];
        }
    };

	typedef IStack<Sample> Samples;

	typedef IIterator<Sample> iterator;
	
	typedef IIterator<const Sample> const_iterator;

    IMPORTANCE_INLINE WeightComparator getWeightComparator(){ return WeightComparator(); }

	IMPORTANCE_INLINE void add(const Importance::Particle & particle, const Vector3 & /*normal*/, const float weight)
	{
		Sample s;
        s.distance = particle.distance;
        IMPORTANCE_ASSERT( isReal( s.distance ) );
		s.original = particle.incidentDir;
		s.photonStart = m_hit.position + particle.incidentDir * particle.distance;
		float dist = (m_hit.position - particle.position).square();
		if ( dist > m_furthest ) m_furthest = dist;				
        Vector2 p;	
#ifdef LIBIMP_DEBUG
        bool res = TDistribution::getMapping().toSquare( m_localFrame, particle.incidentDir, p );
        IMPORTANCE_ASSERT( res );
#else
        TDistribution::getMapping().toSquare( m_localFrame, particle.incidentDir, p );
#endif
		IMPORTANCE_ASSERT( p.x >= 0.0f && p.x <= 1.0f );
		IMPORTANCE_ASSERT( p.y >= 0.0f && p.y <= 1.0f );
		s.shifted_pos = s.pos = QuadVector2( TScalar( p.x ), TScalar( p.y ) );
        //IMPORTANCE_ASSERT(weight > 0.f);
        IMPORTANCE_ASSERT(isReal(weight));
		s.weight = TScalar( weight );
		m_samples.push( s );
	}

    IMPORTANCE_INLINE void clear() {
        m_samples.clear();
    }

	IMPORTANCE_INLINE void reproject( const Vector3 & move )
	{
		for ( Samples::iterator it = m_samples.begin(); it != m_samples.end(); ++it )
		{
#ifdef LIBIMP_DEBUG
            const Vector3 originalDir = (it->photonStart-m_hit.position).getNormalized();
            const Vector3 rotated = m_localFrame.toLocal(originalDir);
            IMPORTANCE_ASSERT(rotated.z >= -0.05f);
#endif

			const Vector3 newPhotonStart = it->photonStart + move;
			const Vector3 wDir = (newPhotonStart - m_hit.position).getNormalized();
			IMPORTANCE_ASSERT( wDir.isReal() && wDir.isNormalized() );
            Vector2 p;			
            TDistribution::getMapping().toSquare( m_localFrame, wDir, p );
			IMPORTANCE_ASSERT( p.x >= 0.0f && p.x <= 1.0f );
			IMPORTANCE_ASSERT( p.y >= 0.0f && p.y <= 1.0f );
			it->shifted_pos = QuadVector2( TScalar( p.x ), TScalar( p.y ) );
		}
	}

	IMPORTANCE_INLINE void prepare( size_t nComponents )
	{
		const int step = (int)nComponents / DIMENSION;
		const int step2 = 2 * step;
		m_pdfs.resize( 2 * step * m_samples.size(), TScalar( 0.0f ) );
		TScalar * iter = &m_pdfs[ 0 ];
		size4 = ( ( m_samples.size() + 3 ) / 4 ) * 4;
		denom = (float * )_aligned_malloc( size4 * sizeof( float ), 16 );
		weight = (float * )_aligned_malloc( size4 * sizeof( float ), 16 );
		size_t i = 0;
		for ( Samples::iterator it = m_samples.begin(); it != m_samples.end(); ++it, iter += step2, ++i )
		{
			it->gamma = iter;
			it->gamma_w = iter + step;
			weight[ i ] = it->weight[ 0 ];
		}
		for (  ; i < size4; ++i )
		{
			denom[ i ] = 1; // log(denom) = log(1) = 0
			weight[ i ] = 0;
		}
	}

	IMPORTANCE_INLINE iterator begin()
	{
		return m_samples.begin();
	}

	IMPORTANCE_INLINE iterator end()
	{
		return m_samples.end();
	}

	IMPORTANCE_INLINE const_iterator begin() const
	{
		return m_samples.cbegin();
	}

	IMPORTANCE_INLINE const_iterator end() const
	{
		return m_samples.cend();
	}

	IMPORTANCE_INLINE size_t size() const
	{
		return m_samples.size();
	}

    IMPORTANCE_INLINE const Sample & operator[]( size_t i ) const 
    {
        IMPORTANCE_ASSERT( i < m_samples.size() );
        return m_samples[ i ];
    }

    IMPORTANCE_INLINE Sample & operator[]( size_t i ) 
    {
        IMPORTANCE_ASSERT( i < m_samples.size() );
        return m_samples[ i ];
    }

	IMPORTANCE_INLINE size_t Size4() const
	{
		return size4;
	}

	IMPORTANCE_INLINE float * Denom() const
	{
		return denom;
	}

	IMPORTANCE_INLINE float * Weight() const
	{
		return weight;
	}

	IMPORTANCE_INLINE float FurthestParticle() const
	{
		return sqrtf( m_furthest );
	}
private:
	Samples m_samples;
	IStack<TScalar> m_pdfs;
	Frame m_localFrame;
	const Hit& m_hit;
	float m_clampDistance;
	float m_furthest;
	float * denom;
	float * weight;
	size_t size4;
};

}

#endif