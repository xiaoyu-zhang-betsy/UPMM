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
#include "Vmf.h"

namespace Importance {    
    namespace EM {

        //////////////////////////////////////////////////////////////////////////
        /// Stepwise E-M adapters
        //////////////////////////////////////////////////////////////////////////

        class VmfmSufficientStats {
            void operator=(const VmfmSufficientStats&);
        public:
            VmfmSufficientStats( VMFM & dist ) : m_dist( dist ) {}

            IMPORTANCE_INLINE void updateWeight( Float /*eta*/, Float /*weight*/ ) {
                //TODO: add weight computing as for gaussians...
#pragma warning(disable:4127)
                IMPORTANCE_ASSERT( false );
                throw std::runtime_error( "progressive E-M for vMFm is not yet implemented" );
            }

            IMPORTANCE_INLINE void updateComponent( int h, Float eta, const Vector3 & dir, Float weight, Float posterior, Float /*invDist*/ ) {       
                //TODO: keep track of harmonic distance
                IMPORTANCE_ASSERT( h < m_dist.nComponents() );

                Vector3 & statsMi   = m_dist.c[ h ].statsMi;
                Float & statsW      = m_dist.c[ h ].statsW;
                
                statsMi = ( 1.f - eta ) * statsMi + eta * weight * posterior * dir;
                statsW  = ( 1.f - eta ) * statsW  + eta * weight * posterior;

                if ( statsW < 1e-7f ) {
                    statsW = 1e-7f;
                }

                IMPORTANCE_ASSERT( !Importance::isNaN( statsW ) && !Importance::isINF( statsW ) );
            }

            IMPORTANCE_INLINE unsigned int getCount() const { return m_dist.k; }

            IMPORTANCE_INLINE void increment() { m_dist.k ++; }
        private:
            VMFM & m_dist;
        };


        class VmfmMixtureModel {
            void operator=(const VmfmMixtureModel&);
        public:
            VmfmMixtureModel( VMFM & dist, Float maxKappa ) : m_maxKappa( maxKappa), m_dist( dist ) {}

            IMPORTANCE_INLINE void update( int h, ERegularizationType /*regularizationType*/ ) {            
                m_dist.c[ h ].computeParams( m_dist.c[ h ].statsMi, m_dist.c[ h ].statsW, m_maxKappa );
                m_dist.w[ h ] = m_dist.c[ h ].statsW;
            }

            IMPORTANCE_INLINE int nComponents() const { return (int) m_dist.nComponents(); }

            IMPORTANCE_INLINE VmfmSufficientStats getStats() { return VmfmSufficientStats( m_dist ); }

            IMPORTANCE_INLINE void setDatasetSize( unsigned int /*val*/ ) {}

            IMPORTANCE_INLINE Float responsibility( int h, const Vector3 & dir ) {
                IMPORTANCE_ASSERT( h < m_dist.nComponents() );
                return m_dist.w[ h ] * m_dist.c[ h ].pdf( dir );
            }

            /// debug
            bool check() const {
                Float error;
                bool isKappasValid;
                Float weightsSum;
                bool isValid = m_dist.check( error, isKappasValid, weightsSum );
                return isValid || !isKappasValid; // don't report kappas above threshold
            }

            std::string toString() const { return m_dist.toString(); }
        private:
            VMFM & m_dist;
            Float m_maxKappa;
        };

    }
}