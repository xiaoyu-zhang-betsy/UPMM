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

#include "Vector3.h"
#include <iostream>


namespace Importance {

    struct Particle {
        Vector3 incidentDir;
        Vector3 position;
        Vector3 normal;
        Float weight;
        Float distance;

        IMPORTANCE_INLINE void serialize( std::ostream & output ) const {
            incidentDir.serialize( output );
            position.serialize( output );
            normal.serialize( output );
            output.write( (char *) &weight, sizeof( Float ) );
            output.write( (char *) &distance, sizeof( Float ) );
        }

        IMPORTANCE_INLINE void deserialize( std::istream & input ) {
            incidentDir.deserialize( input );
            position.deserialize( input );
            normal.deserialize( input );
            input.read( (char *) &weight, sizeof( Float ) );
            input.read( (char *) &distance, sizeof( Float ) );
        }

        IMPORTANCE_INLINE float getWeight(const float distance, const float maxDistance) const {
            const float result = (3/PI)*((maxDistance-distance)/maxDistance);
            IMPORTANCE_ASSERT(distance <= maxDistance);
            return result;
        }
    };
}