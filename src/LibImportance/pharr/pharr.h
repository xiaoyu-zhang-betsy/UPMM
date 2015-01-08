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
#pragma warning(disable:4100)
#pragma warning(disable:4265)

namespace Importance {    

    /************************************************************************/
    /*  Pharr */
    /************************************************************************/
    
    class Pharr : public DefaultDistributionModel {
    public:
        struct Record {
            Vector3 dir;
            Float weight;            
            
            Record() : weight( FLOAT_NAN ) {                
            }

            Record( const Particle & particle ) {
                fromParticle( particle );
            }

            inline void fromParticle( const Particle & particle ) {
                dir     = particle.incidentDir;
                weight  = particle.weight;
            }
        };

    public:
        Pharr() {            
            // this parameter could be static but we don't want to have a cpp file just because of that...
            m_cosAngleDir   = FLOAT_NAN;
            nParticles      = 0;
            m_sum           = 0;            
        }

        void setNormal( const Vector3 & normal ) {
            m_normal = normal;
        }
        
        /*  Distribution interface                                              */
        //////////////////////////////////////////////////////////////////////////

		Float gatherAreaPdfGMM(Vector3 wo, Float radius, Vector2* componentCDFs, Vector2* componentBounds, int &topComponentCDFs, int &topComponentBounds, int baseCDFs, int baseBounds) const{
			IMPORTANCE_ASSERT(false);
			return 1.f;
		}
		Vector3 sampleGatherAreaGMM(Vector2 samples, Vector3 wo, Float radius, int ptrNode, Vector2* componentCDFs, Vector2* componentBounds) const{
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
            for ( int i = 0; i < nParticles; ++i ) {
                Float dotp = Importance::dot( m_data[ i ].dir, dir );
                if ( dotp > 0.999f * m_cosAngleDir ) {
                    pdf += m_data[ i ].weight * squareToUniformConePdf( m_cosAngleDir );				
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

            Vector3 newDirLocal = squareToUniformCone( m_cosAngleDir, Vector2( reusedSample, sample.y ) );		
            return Importance::Frame( m_data[ particleIndex ].dir ).toWorld( newDirLocal );
        }

        //////////////////////////////////////////////////////////////////////////

        inline void add( const Importance::Particle & particle, const Importance::Vector3 & normal, const Float weight ) {
            IMPORTANCE_ASSERT( nParticles < MAX_KNN_PARTICLES );
            Record p( particle );
            m_data[ nParticles++ ] = p;
            m_sum += p.weight;
        }

        inline void setCosAngleDir( Float cosAngleDir ) {
            IMPORTANCE_ASSERT( cosAngleDir >= 0.f && cosAngleDir <= 1.f );
            m_cosAngleDir = cosAngleDir;
        }

        inline int size() const {
            return nParticles;
        } 
        void resize(const int size) {
            this->m_data.resize(size);
        }

		void virtual release() {
			m_data.dealloc();
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