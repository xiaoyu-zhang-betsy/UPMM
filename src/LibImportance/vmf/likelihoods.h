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

#include "datasetrecord.h"
#include "Vmf.h" 
#include "..\shared\Config.h"


namespace Importance {
    class LikelihoodFunctional {
    public:
        LikelihoodFunctional( const VMFM * dist ) : m_dist( dist ) {}

        Float eval( const std::vector<DatasetRecord> & ) const {
            throw std::runtime_error( "not implemented yet" );
            return 0.f;
        }
    protected:
        const VMFM * m_dist;
    };

    class VerbeekLikelihoodFunctional : public LikelihoodFunctional {
    public:
        VerbeekLikelihoodFunctional( const VMFM * dist ) : LikelihoodFunctional( dist ) {};

        Float eval( const Importance::DatasetRecord & rec ) const {
            const VMFM & dist = *m_dist;
            Float denom = 0.f;          
            size_t nComponents = dist.nComponents();
            Float * pst = (Float*) alloca( nComponents * sizeof( Float ) );
            Float res = 0.f;

            for ( size_t h = 0; h < nComponents; ++h ) {
                pst[ h ] = dist.w[ h ] * dist.c[ h ].pdf( rec.dir );
                denom += pst[ h ];
            }

            Float dirWeightDenom = rec.value / denom;
            for ( size_t h = 0; h < nComponents; ++h ) {
                Float posterior = dirWeightDenom * pst[ h ];
                Float contrib = posterior * ( 
                    std::log( dist.w[ h ] / posterior ) + dist.c[ h ].logAveraged( rec.dir )
                    );
                if ( posterior < std::numeric_limits<Float>::min() ) {
                    contrib = 0.f;
                }
                res += contrib;
            }            
            
            return res;
        }

        IMPORTANCE_INLINE Float eval( const std::vector<Importance::DatasetRecord> & dataset ) const {
            Float res = 0.f;                        
            
            for ( size_t i = 0; i < dataset.size(); ++i ) {
                res += eval( dataset[ i ] );
            }

            return res;
        }
    };
}