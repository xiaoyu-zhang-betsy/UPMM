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

#include "Vmf.h"
#include "vmftypes.h"

namespace Importance {   

    template<class TMixtureModel, class TVector>
    void stepwiseEM ( 
        const VmfTypes::Dataset & dataset,
        int nDatasetSize,        
        size_t /*maxIter*/,
        Float alpha,
        TMixtureModel & model, const TVector & ) { 

        IMPORTANCE_ASSERT( alpha > 0.5f && alpha <= 1.f );
                        
        for ( int i = 0; i < nDatasetSize; ++i ) {
            /// compute posterior probabilities
            StaticArray<Float,MAX_FITTED_LOBES> pdf;
            std::fill( &pdf[ 0 ], &pdf[ model.nComponents() ], 0.f );
            Float denom                 = 0.f;;
            const Float & dataWeight    = dataset[ i ].value;            
            const TVector & x           = dataset[ i ].dir;

            for ( int h = 0; h < model.nComponents(); ++h ) {
                pdf[ h ] = dataWeight * model.responsibility( h, x );
                denom += pdf[ h ];
            }
            
            auto stats = model.getStats();
            Float eta = std::pow( stats.getCount() + 1.f, -alpha );
            for ( int h = 0; h < model.nComponents(); ++h ) {
                //////////////////////////////////////////////////////////////////////////
                /// E-step: update sufficient statistics
                //////////////////////////////////////////////////////////////////////////
                Float posterior = pdf[ h ] / denom;
                stats.update( h, eta, x, posterior );

                //////////////////////////////////////////////////////////////////////////
                /// M-step: update model
                //////////////////////////////////////////////////////////////////////////         
                int delayedUpdate = 0;
                if ( stats.getCount() + 1 > delayedUpdate ) {
                    model.update( h );
                }
            }            
            stats.increment();

            IMPORTANCE_ASSERT( model.check() );
        }
    }

    //template<class TMixtureModel, class TVector>
    //void stepwiseEMBatch ( 
    //    const VmfTypes::Dataset & dataset,
    //    int nDatasetSize,        
    //    size_t maxIter,
    //    Float alpha,
    //    TMixtureModel & model, const TVector & ) { 

    //        IMPORTANCE_ASSERT( alpha > 0.5f && alpha <= 1.f );

    //        for ( int i = 0; i < nDatasetSize; ++i ) {
    //            /// compute posterior probabilities
    //            StaticArray<Float,MAX_FITTED_LOBES> pdf;
    //            std::fill( &pdf[ 0 ], &pdf[ model.nComponents() ], 0.f );
    //            Float denom                 = 0.f;;
    //            const Float & dataWeight    = dataset[ i ].value;            
    //            const TVector & x           = dataset[ i ].dir;

    //            for ( int h = 0; h < model.nComponents(); ++h ) {
    //                pdf[ h ] = dataWeight * model.responsibility( h, x );
    //                denom += pdf[ h ];
    //            }
    //        }

    //        auto stats = model.getStats();
    //        Float eta = std::pow( stats.getCount() + 1.f, -alpha );
    //        for ( int h = 0; h < model.nComponents(); ++h ) {
    //            //////////////////////////////////////////////////////////////////////////
    //            /// E-step: update sufficient statistics
    //            //////////////////////////////////////////////////////////////////////////
    //            Float posterior = pdf[ h ] / denom;
    //            stats.update( h, eta, x, posterior );

    //            //////////////////////////////////////////////////////////////////////////
    //            /// M-step: update model
    //            //////////////////////////////////////////////////////////////////////////         
    //            int delayedUpdate = 0;
    //            if ( stats.getCount() + 1 > delayedUpdate ) {
    //                model.update( h );
    //            }
    //        }            
    //        stats.increment();

    //        IMPORTANCE_ASSERT( model.check() );
    //        
    //}
}