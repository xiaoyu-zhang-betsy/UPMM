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

#include <vector>
#include <sstream>

#include "..\LibImportanceTypes.h"
#include "../distribution/DefaultDistributionModel.h"
#include "..\shared\Utils.h"
#include "..\shared\frame.h"
#include "..\viz\AtomicInt.h"

class ImportanceDistribution;

namespace Importance {

    class Jensen : public DefaultDistributionModel {
    private:
        static AtomicInt64 s_nAllocated;
        static size_t s_memory;
    public:
        struct Bin {
            struct Bin() : count( 0 ), val( 0.f ) {}
            int count;
            Float val;
        };
    public:
        Jensen() {
            m_h = m_w = 0;
#ifdef LIBIMP_DEBUG
            used = false;
#endif
        }

#ifdef LIBIMP_DEBUG
        mutable bool used;
#endif

        static void initStatic() {
            s_nAllocated = 0;
            s_memory = 0;
        }

        static inline std::string debugString() {
#ifdef LIBIMP_DEBUG
            std::ostringstream ostr;
            ostr << "Currently allocated: " << s_nAllocated << std::endl
                 <<  "Total memory if not allocated: " << s_memory << std::endl;

            return ostr.str();
#else
            return "Debugging of Jensen's distribution releasing is currently switched off.";
#endif
        }

        void init( size_t w, size_t h, const Vector3 & normal ) {
            m_hist.resize( w * h );        
            m_normal                = normal;
            m_w                     = w;
            m_h                     = h;
            m_totalWeight           = 0.f;
            m_isInitialized         = true;
            m_isPreparedForSampling = false;
            m_count                 = 0;

#ifdef LIBIMP_DEBUG
            ++Jensen::s_nAllocated;            
            Jensen::s_memory += sizeof( m_hist ) + w * h * sizeof( *this );
#endif
        }

        bool prepareForSampling() {
            if ( size() == 0 ) {
                return m_isPreparedForSampling = false; 
            }

            m_rowsSum.clear();
            m_rows.clear();

            m_rowsSum.reserve( m_h );
            m_rows.resize( m_h );            
            
            for ( size_t j = 0; j < m_h; ++j ) {
                Float rowSum = 0.f;
                m_rows[ j ] = DiscreteDistribution( m_w );
                for ( size_t i = 0; i < m_w; ++i ) {
                    Float val = m_hist[ getIndex( i, j ) ].val;
                    rowSum += val;
                    m_rows[ j ].append( val );
                }
                m_rows[ j ].normalize();
                m_rowsSum.append( rowSum );
            }

            m_rowsSum.normalize();

            IMPORTANCE_ASSERT( internalCheck() );

            s_memory += m_h + sizeof( Float ) + m_h * sizeof( DiscreteDistribution ) + m_w * m_h * sizeof( Float );

            return m_isPreparedForSampling = true;
        }

        bool internalCheck() {
            if ( size() == 0 ) {
                return true;
            }

            Float sum = 0.f;
            for ( size_t j = 0; j < m_h; ++j ) {
                for ( size_t i = 0; i < m_w; ++i ) {
                    sum += m_hist[ getIndex( i, j ) ].val;
                }
            }

            if ( sum == 0.f ) {
                return false;
            }

            Float diff = std::abs( 1.f - m_totalWeight / sum );
            return diff < 1e-5;
        }

        Float pdf( const Vector3 & wDir ) const {
#ifdef LIBIMP_DEBUG
            used = true;
#endif
            IMPORTANCE_ASSERT(m_isInitialized);
            Vector3 lDir = Frame( m_normal ).toLocal( wDir );            
        
            if ( lDir.z < 0.f ) { 
                return 0.f; 
            }

            if ( size() == 0 ) {                
                Float res = squareToCosineHemispherePdf( lDir );
                IMPORTANCE_ASSERT( res >= 0.f );
                return res;                
            }

            IMPORTANCE_ASSERT( m_isInitialized );
            IMPORTANCE_ASSERT( m_totalWeight > 0.f );
            
            Float uvPdf = m_hist[ getIndex( lDir ) ].val / m_totalWeight * m_w * m_h;            
            //Float sinTheta = std::sqrt( std::max( 0.f, 1.f - lDir.z * lDir.z ) );
            // Jacobian of Shirley mapping from [0,1]^2 to hemisphere is 2*pi
            Float jacobi = 2.0f * IMP_PI; // sinTheta * IMP_PI * IMP_PI;            
            Float res = uvPdf / jacobi;
            IMPORTANCE_ASSERT( res >= 0.f );
            return res;
        }

		Float gatherAreaPdf(Vector3 wo, Float radius, std::vector<Vector2> &componentCDFs, std::vector<Vector2> &componentBounds, int baseCDFs, int baseBounds) const{
			IMPORTANCE_ASSERT(false);
			return 1.f;
		}

        Importance::Vector3 sampleDirection( const Importance::Vector2 & sample ) const {
#ifdef LIBIMP_DEBUG
            used = true;
#endif
            if ( size() == 0 ) {
                Vector3 lDir = squareToCosineHemisphere( sample );
                IMPORTANCE_ASSERT( lDir.z > 0.f );
                return Frame( m_normal ).toWorld( lDir );
            }

            IMPORTANCE_ASSERT( m_isInitialized );
            IMPORTANCE_ASSERT( m_isPreparedForSampling );

            Vector2 tmp( sample );
        
            size_t j = m_rowsSum.sampleReuse( tmp.y ),
                   i = m_rows[ j ].sampleReuse( tmp.x );

            Float u = ( i + tmp.x ) / m_w,
                  v = ( j + tmp.y ) / m_h;

            IMPORTANCE_ASSERT( u >= 0.f && u <= 1.f && v >= 0.f && v <= 1.f );
            Vector3 lDir = uniformHemisphereShirley( u, v );
            //Vector3 lDir = squareToHemisphere( u, v );
            return Frame( m_normal ).toWorld( lDir );            
        }

        virtual void release() {
            IMPORTANCE_ASSERT(used);
            for ( auto it = m_rows.begin(); it < m_rows.end(); ++it ) {
                it->release();
            }
            std::vector<DiscreteDistribution>().swap( m_rows );
            std::vector<Bin>().swap( m_hist );
            m_rowsSum.release();
            m_isPreparedForSampling = false;
            m_isInitialized = false;
#ifdef LIBIMP_DEBUG
            --Jensen::s_nAllocated;
#endif
        }

        virtual std::string toString() const {
            std::ostringstream ostr;
            ostr << "Jensen [ " << std::endl
                 << "w x h = " << m_w << " x " << m_h << std::endl
                 << "]" << std::endl;
            return ostr.str();
        }

        //////////////////////////////////////////////////////////////////////////

        inline void add( const Importance::Particle & particle, const Importance::Vector3 & /*normal*/, const Float /*weight*/ ) {        
            IMPORTANCE_ASSERT( m_isInitialized );
            IMPORTANCE_ASSERT( particle.weight > 0.f );
            Vector3 lDir = Frame( m_normal ).toLocal( particle.incidentDir );
            IMPORTANCE_ASSERT( lDir.z >= 0.f );                        
            m_hist[ getIndex( lDir ) ].val      += particle.weight;
            m_hist[ getIndex( lDir ) ].count ++;
            m_totalWeight                       += particle.weight;
            m_isPreparedForSampling             = false;
            m_count ++;
        }

        inline size_t size() const { return m_count; }

    protected:
        size_t getIndex( const Vector3 & lDir ) const {
            IMPORTANCE_ASSERT( m_w > 0 && m_h > 0 );
            Float u, v;            
            uniformHemisphereShirleyInverse( lDir, u, v );
            //hemisphereToSquare( lDir, u, v );
            size_t x = std::min( (size_t) std::floor( u * m_w ), m_w - 1 ),
                   y = std::min( (size_t) std::floor( v * m_h ), m_h - 1 );
            return getIndex( x, y );
        }

        size_t getIndex( size_t x, size_t y ) const {
            IMPORTANCE_ASSERT( x >= 0 && x < m_w && y >= 0 && y < m_h );
            return y * m_w + x;
        }

    private:
        std::vector<Bin> m_hist;
        size_t m_w, m_h;
        Vector3 m_normal;
        Float m_totalWeight;
        DiscreteDistribution m_rowsSum;
        std::vector<DiscreteDistribution> m_rows;
        bool m_isPreparedForSampling, m_isInitialized;
        size_t m_count;
        friend class ::ImportanceDistribution;
    };

}