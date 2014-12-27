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

#include <fstream>

namespace Importance {
    inline void dumpNNParticles( const Hit & hit, const ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> & nnqueryRes, const IStack<Particle> & particles, int found ) {        
        if ( hit.position.isReal() ) {
            std::ofstream file( "particles.dump", std::ofstream::trunc | std::ofstream::binary );
            hit.serialize( file );
                        
            for ( int i = 0; i < found; i++ ) {
                const Particle & part = particles[ nnqueryRes[ i ].index ];
                part.serialize( file );
            }

            file.close();
        }
    }    

    inline int loadNNParticles( const std::string & filename, Hit & hit, ImStaticArray<Particle, MAX_KNN_PARTICLES> & particles ) {
        int nParticles = 0;
        std::ifstream file;
        file.exceptions ( std::ifstream::failbit | std::ifstream::badbit );        
        try {
            file.open( filename, std::istream::binary );
            hit.deserialize( file );

            while ( file.peek() != EOF && nParticles < MAX_KNN_PARTICLES ) {
                particles[ nParticles++ ].deserialize( file );
            }
        } catch ( std::ifstream::failure e ) {
            ILog( EWarn, "Error when opening the file " + filename + " - " + std::string( e.what() ) );
            return 0;
        }

        return nParticles;
    }
}