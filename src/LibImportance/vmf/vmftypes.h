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
    namespace VmfTypes {
        typedef ImStaticArray<DatasetRecord,MAX_KNN_PARTICLES> Dataset;

        class DatasetAdapter {
            void operator=(const DatasetAdapter&);
        public:
            DatasetAdapter( Dataset & dataset ) 
                : m_dataset( dataset ), m_particlesUsed( 0 ), m_totalWeight( 0.f ) {}

            IMPORTANCE_INLINE void add( const Particle & particle, const Vector3 & /*normal*/, const Float weight ) {
                m_dataset[ m_particlesUsed ].dir       = particle.incidentDir;
                m_dataset[ m_particlesUsed ].invDist   = 1.f/particle.distance;
                m_dataset[ m_particlesUsed ].value     = particle.weight;
                m_totalWeight += weight;
                ++m_particlesUsed;                
            }
            int    getNParticlesUsed() const { return m_particlesUsed; }
            Float  getTotalWeight() const { return m_totalWeight; }
        private:
            Dataset & m_dataset;
            int m_particlesUsed;
            Float m_totalWeight;
        };
    }
}