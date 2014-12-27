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


#include "Vmf.h"
#include "vmfmfactory.h"
#include "..\shared\Utils.h"
#include "..\shared\FastFloat\FastFloat.h"
#include "..\shared\simplelogger.h"
#include <iostream>

#pragma warning(push)
#pragma warning(disable:4127)
#pragma warning(disable:4100)

namespace Importance {

    void VMFM::randomInitHemi( Importance::Random & rnd, const Vector3 & nrm ) {
        Frame frame( nrm );
        for ( size_t h = 0; h < nComponents(); ++h ) {        
            Vector3 localDir = squareToHemisphere( Importance::nextFloat(rnd), Importance::nextFloat(rnd) );
            this->c[ h ].m_mi = frame.toWorld( localDir );
            this->c[ h ].m_kappa = 1.0f;
            this->w[ h ] = 1.f / nComponents();
            //we suppose that this method is used during EM algorithm
            this->c[ h ].initPDFOnly();
        }
        //cdf is still not computed - we might not need that (if we don't aim to sample from this dist)
        m_initialized = false;
    }

    size_t VMFM::getSampledComponentIndex( Float u ) const {
        IMPORTANCE_ASSERT( u >= 0.f && u <= 1.f );
        //inverse cdf function            
        size_t i = 0;
        for( ; u >= cdf[ i ] && i < nComponents() + 1; ++i );
        return std::min( i - 1, (size_t)nComponents() - 1 );
    }


    const Float VMF::UNIFORM_PDF = 1.0f / (4.0f * IMP_PI);
    const Float VMF::TWOPI = 2.0f * IMP_PI;
    const Float VMF::EXP_MINUSTWO = std::exp( -2.0f );

    const Float VMFM::EPSILON   = 1e-3f;
    const Float VMFM::MAX_KAPA  = 5e6f;    
    const Float VMFM::MINIMAL_ERROR = 1e-5f;

}

const Importance::Float Importance::VectorFloatTraits::initVal = Importance::NAN;


#pragma warning(pop)