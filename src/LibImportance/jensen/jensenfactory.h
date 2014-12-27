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

#include "../shared/basicfactory.h"
#include "jensen.h"

namespace Importance {

    class JensenFactory : public BasicDistributionFactory {
    public:
        virtual bool getInstance( 
            const Importance::Hit& hit, 
            const Importance::KdQueryResult* particles,
            int nParticlesFound, 
            DefaultDistributionModel& output ) {
                IMPORTANCE_ASSERT( &dynamic_cast<Jensen&>( output ) );
                Jensen & dist = (Jensen& ) output;
                dist.init( m_config->jensen.w, m_config->jensen.h, hit.normal );
                DataFilter<Jensen> particleFilter( m_photonMap, *m_config );
                particleFilter.prepareData( particles, hit, nParticlesFound, dist.m_cacheStats, dist );
                dist.prepareForSampling();
                return true;
        }


        virtual bool getInstanceProgressive( 
            const Hit& hit, 
            const KdQueryResult* particles,
            int nParticlesFound, 
            DefaultDistributionModel & output ) {
                IMPORTANCE_ASSERT( &dynamic_cast<Jensen&>( output ) );                
                Jensen & dist = (Jensen& ) output;                
                DataFilter<Jensen> particleFilter( m_photonMap, *m_config );
                particleFilter.prepareData( particles, hit, nParticlesFound, dist.m_cacheStats, dist );
                dist.prepareForSampling();
                return true;
        }

        virtual std::string toString() const {
            std::ostringstream strs;
            strs << "JensenFactory [ " << std::endl                               
                << "w = " << m_config->jensen.w << std::endl
                << "h = " << m_config->jensen.h << std::endl
                << "]" << std::endl;
            return strs.str();
        }
    };
}