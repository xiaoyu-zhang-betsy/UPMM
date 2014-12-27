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
#include "../shared/Vector3.h"
#include "../shared/StaticArray.h"
#include "../shared/Stack.h"
#include "../shared/PointKdTreeShared.h"

#pragma warning(push)
#pragma warning( disable : 4512 )

namespace Importance {
    class IDensityEstimator {
    public:
        virtual Float computeDensity( const ImStaticArray<KdQueryResult, ENVIRO_DIR_SEARCH_COUNT> & nearest, int found, Float sqrRadius ) const = 0;
        virtual void knnQuery( const Vector3 & position, int minParticles, int & found, Float & sqrSearchRad, 
            ImStaticArray<KdQueryResult, ENVIRO_DIR_SEARCH_COUNT> & nearest ) const = 0;
    };

    class PDensityEst {
    private:
        IStack<Float> m_reduction;
        int m_depth;
        IDensityEstimator & m_densityEstimator;
        int m_countUpdates;
        int m_nIter;        
    public:
        struct Pixel {
            //Float sqrRefRad;
            Float value;
            Vector3 position;            
#ifdef LIBIMP_DEBUG
            Pixel() : value( NAN ), position( NAN )/*, sqrRefRad( NAN )*/ {}
#endif
        };

        PDensityEst( Float alpha, IDensityEstimator & densityEstimator, int depth = 80 ) 
            : m_densityEstimator( densityEstimator ), m_depth ( depth ), m_countUpdates( 0 ), m_nIter( 0 ) {
            init( alpha );
        }      

        IMPORTANCE_INLINE int getCountUpdates() const {
            return m_countUpdates;
        }

        IMPORTANCE_INLINE void init( Float alpha ) {
            m_reduction.resize( m_depth );
            m_reduction[ 0 ] = 1.f;
            for ( int i = 1; i < m_depth; i++ ) {
                m_reduction[ i ] = m_reduction[ i - 1 ] * ( i + alpha ) / (Float) ( i + 1 );
            }   
        }

        IMPORTANCE_INLINE void estimateKnnRefRad( Pixel & pixel ) {
            const int & i = m_nIter;            
            IMPORTANCE_ASSERT( i < m_depth );
            IMPORTANCE_ASSERT( isReal( pixel.value ) );
            if ( i >= m_depth ) {
                throw std::runtime_error( "Basic depth of precomputed radius reduction was exceeded" );
                return;
            }            

            ImStaticArray<KdQueryResult, ENVIRO_DIR_SEARCH_COUNT> nearest;
            Float sqrRefRad;
            int found;
            const int minParticles = 15;
            m_densityEstimator.knnQuery( pixel.position, minParticles, found, sqrRefRad, nearest );
            Float sqrRad = m_reduction[ i ] * sqrRefRad;                        
            
            // if there is less than minimum particles, we assume kernel density query
            Float furthestSqrRad = sqrRefRad;
            if ( found >= minParticles ) {
                // We have enough particles to shrink the radius
                filter( sqrRad, minParticles, nearest, found );
                furthestSqrRad = 0.f;
                for ( int i = 0; i < found; ++i ) {
                    furthestSqrRad = std::max( nearest[ i ].distSqr, furthestSqrRad );
                }
            }

            pixel.value = pixel.value * i + m_densityEstimator.computeDensity( nearest, found, furthestSqrRad );
            //only knn query
            //pixel.value = m_densityEstimator.computeDensity( nearest, found, sqrRad );
            pixel.value /= ( i + 1 );
            IMPORTANCE_ASSERT( isReal( pixel.value ) );
            ++m_countUpdates;
        }

        IMPORTANCE_INLINE int nextIteration() {
            m_countUpdates = 0;
            return ++m_nIter;
        }

        protected:
            /* Get rid off particles that are further than maxSqrRad but preserve at least minCount of them. 
               The furthest particles goes first. 
               !!!Important!!!
               Results input must have a max-heap property.
               */
            IMPORTANCE_INLINE void filter( const Float maxSqrRad, int minCount, ImStaticArray<KdQueryResult, ENVIRO_DIR_SEARCH_COUNT> & results, int & size ) {
                if ( size <= minCount ) {
                    return;
                }

                ImExternalHeap<KdQueryResult, KdQueryResultComparator> heap( results.ptr(), ENVIRO_DIR_SEARCH_COUNT, size );
                while( heap.peek().distSqr > maxSqrRad && heap.size() > minCount ) {
                    heap.removeTop();
                }
                size = heap.size();
            }
    };
}

#pragma warning(pop)