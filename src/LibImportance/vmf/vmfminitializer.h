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
#include "../LibImportanceTypes.h"
#include "random.h"
#include "datasetrecord.h"
#include "Vmf.h"
#include "greedy.h"
#include "vmftypes.h"

#pragma warning(push)
#pragma warning(disable:4100)

namespace Importance {

    class VmfmInitializer {
    public:
        IMPORTANCE_INLINE void init(
            Config::VMFM::InitType method,
            const VmfTypes::Dataset & dataset,
            size_t nDataSize,            
            const Vector3 & n,
            VMFM & vmfm ) {

            switch ( method ) {                
            case Importance::Config::VMFM::EKMeansPP:                    
                kmeanspp( vmfm, dataset, nDataSize, m_rnd );
                break;

            case Importance::Config::VMFM::ERandom:                    
                vmfm.randomInitHemi( m_rnd, n );        
                break;

            case Importance::Config::VMFM::EPhotonUniform:
                {                    
                    Float initialWeight = 1.f / vmfm.nComponents();
                    for ( size_t h = 0; h < vmfm.nComponents(); ++h ) {            
                        /// simple stratification (no direction should be chosen more than once)
                        Float u = ( h + Importance::nextFloat(m_rnd) ) / ( Float ) vmfm.nComponents();
                        int index           = (int) std::min( nDataSize - 1, (size_t) ( u * nDataSize ) );
                        vmfm.w[ h ]         = initialWeight;
                        vmfm.c[ h ].m_kappa  = 1.0f;
                        vmfm.c[ h ].m_mi    = dataset[ index ].dir;
                        vmfm.c[ h ].initPDFOnly();
                    }
                }
                break;

            case Importance::Config::VMFM::EPhotonWeighed:
                {
                    std::vector<Float> cdf( nDataSize + 1, 0.f );
                    //compute cdf        
                    Float sum = 0.0f;
                    /// cdf is postponed by one index for comfy samples reusing
                    cdf[ 0 ] = 0.f;
                    for( int j = 1; j < nDataSize + 1; ++j ) {            
                        sum += dataset[ j - 1 ].value;
                        cdf[ j ] = sum;
                    }                        

                    Float initialWeight = 1.f / vmfm.nComponents();
                    for ( size_t h = 0; h < vmfm.nComponents(); ++h ) {
                        int i = 0;
                        Float u = Importance::nextFloat(m_rnd) * sum;
                        for( ; u >= cdf[ i ] && i < nDataSize + 1; ++i );
                        int index           = std::min( i - 1, (int) nDataSize - 1 );
                        vmfm.w[ h ]         = initialWeight;
                        vmfm.c[ h ].m_kappa  = 1.0f;
                        vmfm.c[ h ].m_mi    = dataset[ index ].dir;
                        vmfm.c[ h ].initPDFOnly();
                    }

                }
                break;

            case Importance::Config::VMFM::ESingular: 
                {
                    Float initialWeight = 1.f / vmfm.nComponents();
                    for ( size_t h = 0; h < vmfm.nComponents(); ++h ) {
                        vmfm.w[ h ]         = initialWeight;
                        vmfm.c[ h ].m_kappa  = 1.0f;
                        vmfm.c[ h ].m_mi    = n;
                        vmfm.c[ h ].initPDFOnly();
                    }
                }
                break;

            default:
                throw std::runtime_error( "Unexpected case" );
            }                
        }

    private:
        void kmeanspp_version_from_greedy( VMFM & vmfm, const VmfTypes::Dataset & dataset, size_t nDataSize, Importance::Random & rnd ) {
            GreedyToolSet::kmeanspp( &dataset[0], nDataSize, rnd, &vmfm.c[0], &vmfm.w[0], vmfm.nComponents() );
        }

        void kmeanspp( VMFM & vmfm, const VmfTypes::Dataset & dataset, size_t nDataSize, Importance::Random & rnd ) {
                size_t nComponents  = vmfm.nComponents();
                Float initialWeight = 1.f / nComponents;

                /// TODO: optimize this code; generate samples in batches

                /// Chose the first component direction at random
                Float u = Importance::nextFloat(rnd);
                int indx = (int) (nDataSize * u);
                IMPORTANCE_ASSERT( indx < nDataSize );
                vmfm.c[ 0 ].m_mi    = dataset[ indx ].dir;
                vmfm.c[ 0 ].m_kappa  = 1.0f;    
                vmfm.w[ 0 ]         = initialWeight;
                vmfm.m_nComponents  = 1;
                IMPORTANCE_ASSERT( vmfm.c[ 0 ].m_mi.isNormalized() );
                vmfm.c[ 0 ].initPDFOnly();

                Float * pdf  = (Float*) alloca( ( nDataSize + 1 ) * sizeof( Float ) );
                for ( size_t iComp = 1; iComp < nComponents; ++iComp ) {

                    Float sum = 0.f;                                    
                    for ( int i = 0; i < nDataSize; ++i ) {
                        Float maxPdf =  0.f; /// minimum distance
                        for ( size_t j = 0; j < iComp; ++j ) {
                            /// probably this metric doesn't have to be pdf from lobes, could be just a cosine of angles for example
                            /// but I think the exponential shape does match better our problem - there is a bigger penalty for those 
                            /// directions which are too close and won't be chosen that often
                            /// although as an optimization the normalization factor could be removed as it will do no difference
                            maxPdf = std::max( vmfm.c[ j ].pdf( dataset[ i ].dir ), maxPdf );                    
                        }
                        pdf[ i ] = maxPdf; //vmfm.pdf( dataset[ i ].dir );

                        /// we need to know this so that we can easily twist the process of choice from CDF                
                        /// in our metric the highest value (distance) must has the lowest chance of being chosen
                        pdf[ i ]    = 1.f / pdf[ i ];
                        sum        += pdf[ i ];
                    }

                    u = Importance::nextFloat(rnd);
                    u = u * sum;
                    Float cdfVal = pdf[ 0 ];
                    int l;
                    pdf[ nDataSize ] = 0.f;
                    for ( l = 0; l < nDataSize && cdfVal <= u; ++l ) {
                        cdfVal += pdf[ l + 1 ];
                    }
                    IMPORTANCE_ASSERT( l < nDataSize && cdfVal > u );

                    vmfm.c[ iComp ].m_mi    = dataset[ l ].dir;
                    vmfm.c[ iComp ].m_kappa  = 1.0f; 
                    vmfm.w[ iComp ]         = initialWeight;
                    vmfm.m_nComponents     += 1;
                    IMPORTANCE_ASSERT( vmfm.c[ iComp ].m_mi.isNormalized() );
                    vmfm.c[ iComp ].initPDFOnly();
                }
            }
    
    private:
        Random m_rnd;
    };
}

#pragma warning(pop)