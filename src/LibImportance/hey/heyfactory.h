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


/*
    This file contains implementation of importance sampling technique according to article
    Hey, H., Purgathofer, W.: "Importance Sampling with Hemispherical Particle Footprints"

    Author: Jiri Vorba
*/

#pragma once

#include "../shared/basicfactory.h"
#include "hey.h"

#pragma warning(push)
#pragma warning(disable:4100)


namespace Importance {
    namespace HeyFactoryConsts {
        enum { RES_X = 64, RES_Y = 64 };
        const Float kr = 240.f;
    }

    /************************************************************************/
    /* Splatter */
    /************************************************************************/

    /** Computes radius of particles based on their directional density. Uses splatting  
     *  of incident particles to a discrete grid.
     */
    class Splatter {
    public:
        enum { RES_X = HeyFactoryConsts::RES_X, RES_Y = HeyFactoryConsts::RES_Y };
    public:
        Splatter( const Vector3 & normal, const Float (* angles)[ RES_Y ] ) : m_nParticles( 0 ), m_angles( angles ) {
            memset( m_grid, 0, sizeof( int ) * RES_X * RES_Y );
            m_coords = Frame( normal );
        }

        inline void add( const Importance::Particle & particle, const Importance::Vector3 & normal, const Float ) {
            IMPORTANCE_ASSERT( m_nParticles < MAX_KNN_PARTICLES );                                    
            m_particles[ m_nParticles++ ] = particle;
            splat( particle.incidentDir );
        }

        inline void fill( Hey & heyDist ) const {
            for ( int i = 0; i < m_nParticles; ++i ) {        
                Vector3 locaDir = m_coords.toLocal( m_particles[ i ].incidentDir );
                heyDist.add( m_particles[ i ], computeRadius( locaDir ), NAN );
            }
        }

    private:
        inline Float computeRadius( const Vector3 & localDir ) const {
            int x, y;
            toCell( localDir, x, y );

            Float rmax = localDir.z;
            if ( m_angles[ x ][ y ] < 0.f ) {
                IMPORTANCE_ASSERT( rmax > 0.f );
                return rmax;
            }

            IMPORTANCE_ASSERT( m_grid[ x ][ y ] > 0 && m_grid[ x ][ y ] <= m_nParticles );

            Float r = HeyFactoryConsts::kr * m_angles[ x ][ y ] / m_grid[ x ][ y ];
            IMPORTANCE_ASSERT( r > 0.f && rmax > 0.f );
            return std::min( r, rmax );
        }

        inline void splat( const Vector3 & direction ) {
            int sx, sy;
            toCell( m_coords.toLocal( direction ), sx, sy );
            int startX  = std::max( 0, sx - 1 ), 
                endX    = std::min( sx + 1, RES_X - 1 );
            int startY  = std::max( 0, sy - 1 ), 
                endY    = std::min( sy + 1, RES_Y - 1 );
            
            for ( int y = startY; y <= endY; ++y )
                for ( int x = startX; x <= endX; ++x )
                    m_grid[ x ][ y ]++;
        }

        inline void toCell( const Vector3 & localDir, int & sx, int & sy ) const {            
            sx      = (int) std::floor( ( localDir.x + 1.f ) * 0.5f * RES_X );
            sy      = (int) std::floor( ( localDir.y + 1.f ) * 0.5f * RES_Y );

            IMPORTANCE_ASSERT( sx >= 0 && sx < RES_X && sy >= 0 && sy < RES_Y );
        }

    private:
        Importance::ImStaticArray<Particle, MAX_KNN_PARTICLES> m_particles;
        int m_nParticles;      
        Frame m_coords;
        int m_grid[ RES_X ][ RES_Y ];
        const Float (* m_angles)[ RES_Y ];
    };

    /************************************************************************/
    /* HeyFactory     */
    /************************************************************************/

    class HeyFactory : public BasicDistributionFactory { 
    public:
        enum { RES_X = HeyFactoryConsts::RES_X, RES_Y = HeyFactoryConsts::RES_Y };
    public:
        inline void init ( const IStack<Importance::Particle>& photons, const Importance::Config * config ) {
            BasicDistributionFactory::init( photons, config );
            precomputeSolidAngles();
        }

       virtual bool getInstance( 
            const Importance::Hit & hit, 
            const Importance::KdQueryResult * photons,
            int found, 
            DefaultDistributionModel & output ) {
                IMPORTANCE_ASSERT( &dynamic_cast<Hey&>( output ) );
                Hey & heyDist = (Hey& ) output;
                heyDist.setNormal( hit.normal );
                heyDist.resize(found);                
                
                Splatter splatter( hit.normal, m_angles );
                DataFilter<Splatter> particleFilter( m_photonMap, *m_config );
                
                particleFilter.prepareData( photons, hit, found, output.m_cacheStats, splatter );
                splatter.fill( heyDist );               

                return true;
        }

        virtual std::string toString() const {
            std::ostringstream strs;
            strs << "PharrFactory [ " << std::endl                               
                 << "]" << std::endl;
            return strs.str();
        }

    private:
        void precomputeSolidAngles() {            
            for ( int y = 0; y < RES_Y; ++y ) {
                for ( int x = 0; x < RES_X; ++x ) {
                    Float rx = 2.f * ( x + 0.5f ) / RES_X - 1.f,
                          ry = 2.f * ( y + 0.5f ) / RES_Y - 1.f;
                    Float c = safe_sqrt( 1.f - rx * rx - ry * ry );
                    if ( c > 0.f ) {
                        m_angles[ x ][ y ] = 4.f / ( c * RES_X * RES_Y );
                    } else {
                        // indicate discretization problem so that the radius is set to rmax
                        m_angles[ x ][ y ] = -1.f;
                    }
                }
            }
        }
    private:
        Float m_angles[ RES_X ][ RES_Y ]; 
    };

}

#pragma warning(pop)