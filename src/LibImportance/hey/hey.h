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

#include "../LibImportanceTypes.h"
#include "../distribution/DefaultDistributionModel.h"
#include "../shared/Utils.h"
#include "../shared/frame.h"
#include "../shared/StaticArray.h"
#include "../shared/simplelogger.h"

#pragma warning(push)
#pragma warning(disable:4127)

namespace Importance {    

    /************************************************************************/
    /*  Pharr */
    /************************************************************************/

    class Hey : public DefaultDistributionModel {
    public:
        struct Record {
            Vector3 dir;
            Float weight;
            Float r;

            Record() : weight( FLOAT_NAN ), r( FLOAT_NAN ) {                
            }

            Record( const Particle & particle, Float radius ) {
                r = radius;
                fromParticle( particle );
            }

            inline void fromParticle( const Particle & particle ) {
                dir     = particle.incidentDir;
                weight  = particle.weight;
            }
        };

        //typedef Importance::ImStaticArray<Record, Importance::MAX_KNN_PARTICLES> DataArray;
    public:
        Hey() {                                    
            nParticles      = 0;
            m_sum           = 0;
        }

        void setNormal( const Vector3 & normal ) {
            m_normal = normal;
        }

        /*  Distribution interface                                              */
        //////////////////////////////////////////////////////////////////////////

		Float gatherAreaPdf(Vector3 wo, Float radius, std::vector<Vector2> &componentCDFs, std::vector<Vector2> &componentBounds, int baseCDFs, int baseBounds) const{
			IMPORTANCE_ASSERT(false);
			return 1.f;
		}
		Vector3 sampleGatherArea(Vector2 samples, Vector3 wo, Float radius, int ptrNode, std::vector<Vector2> componentCDFs, std::vector<Vector2> componentBounds) const{
			IMPORTANCE_ASSERT(false);
			return Vector3(0.f);
		}

        Float pdf( const Vector3 & dir ) const {                        
            if ( nParticles == 0 ) {
                Vector3 lDir = Frame( m_normal ).toLocal( dir );
                if ( lDir.z < 0.f ) {
                    return 0.f;
                }
                return squareToCosineHemispherePdf( lDir );
            }

            Float pdf = 0.f;
            Float cosAngle;
            for ( int i = 0; i < nParticles; ++i ) {
                Float d = distance( m_data[ i ].dir, dir, cosAngle );
                if ( d < m_data[ i ].r && cosAngle > 0.f ) {                            
                    pdf += m_data[ i ].weight / footprintSolidAngle( m_data[ i ], dir );
                }
            }

            return pdf / m_sum;
        }


        Importance::Vector3 sampleDirection( const Importance::Vector2 & sample ) const {
            if ( nParticles == 0 ) {
                Vector3 lDir = squareToCosineHemisphere( sample );
                IMPORTANCE_ASSERT( lDir.z > 0.f );
                return Frame( m_normal ).toWorld( lDir );
            }

            Float cdfVal = m_data[ 0 ].weight;
            Float u = sample.x * m_sum;
            int i;
            for ( i = 0; ( i < nParticles - 1 ) && ( cdfVal <= u ); ++i ) {
                cdfVal += m_data[ i + 1 ].weight;
            }
            int particleIndex = std::min<int>( i, (int)( nParticles - 1 ) );             

            Float reusedSample = ( cdfVal - u ) / ( m_data[ particleIndex ].weight );
            IMPORTANCE_ASSERT( reusedSample >= 0.f - 1e4 && reusedSample <= 1.f + 1e4 );
            clamp( reusedSample, 0.f, 1.f );
            
            const Record & rec = m_data[ particleIndex ];
            Vector2 uv = squareToUniformDiskConcentric( Vector2( reusedSample, sample.y ) );            
            uv.x *= rec.r; uv.y *= rec.r;
            Vector3 newDirLocal = Vector3( uv.x, uv.y, safe_sqrt( 1.f - uv.x * uv.x - uv.y * uv.y ) );
            return Importance::Frame( rec.dir ).toWorld( newDirLocal );
        }

        //////////////////////////////////////////////////////////////////////////

        inline void add( const Importance::Particle & particle, Float radius, const Float ) {
            IMPORTANCE_ASSERT( nParticles < MAX_KNN_PARTICLES );  
            IMPORTANCE_ASSERT( radius > 0.f );
            Record p( particle, radius );                       
            m_data[ nParticles++ ] = p;
            m_sum += p.weight;
        }

        inline int size() const {
            return nParticles;
        }

        void resize(const int size) {
            m_data.resize(size);
        }

		virtual void release() {
			m_data.dealloc();
		}

    private:
        inline Float distance( const Vector3 & dirp, const Vector3 & dir, Float & cosAngle ) const {
            cosAngle = dot( dir, dirp );
            return ( cosAngle * dirp - dir ).length();
        }

        inline Float footprintSolidAngle( const Record & rec, const Vector3 & ) const {
            IMPORTANCE_ASSERT( rec.r >= 0.f );
            Float rsqr  = rec.r * rec.r;
            Float c     = 1 - safe_sqrt( 1 - rsqr );
            if ( c == 0.f ) {
                return IMP_PI * rsqr;
            }
            
            return 2.f * IMP_PI * c;
        }
    private:
        IStack<Record> m_data;        
        int nParticles;
        Float m_sum;
        Float m_cosAngleDir;
        Vector3 m_normal;
    };
}

#pragma warning(pop)