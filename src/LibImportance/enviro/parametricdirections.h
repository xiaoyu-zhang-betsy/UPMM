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
#include "ienviromap.h"

/* We could avoid this include by forward declaration */
#include "../shared/PointKdTree.h"
#include "../caching/CacheStats.h"

#include "../gaussian/arraysse.h"
#include "../gaussian/gaussiansse.h"
#include "../gaussian/gaussian_implsse.h"
#include "../em/stepwise_emsse.h"
#include "../gaussian/gaussian_stepwise_emsse.h"

namespace Importance {
    class ParametricDirections {
    private:
        static const int DIMENSION      = 4;
        /* Because there are usually tens of thousands of particles on the environment we 
          split them and feed the fitting in smaller batches (thus we don't do dataset adjustment 
          on such a big datasets) 
          For first of those batches, we run the static stepwise EM.
          */
        static const size_t BATCH_SIZE  = 1000;

    public:
        typedef ArrayGaussianSSE< Float4, DIMENSION, GaussianMixtureSSESphere> Array;
        typedef EM::GaussianMixtureModel< Float4, DIMENSION, GaussianMixtureSSESphere> GaussianMixtureModel;

    public:
        ParametricDirections( const IStack<Particle>& particles, 
            const PointKdTree<Kd3PositionTraits> & /*tree*/ )
        : m_particles( particles ) {
            m_dummyHit.normal   = Vector3( 0.f, 0.f, 1.f );
            m_dummyHit.position = Vector3( 0.f );
        }


        IMPORTANCE_INLINE void init( const EnviroMap * enviro, 
            ImBitmap<Float> & raster, DebuggingBundle & debug, const Config & cfg ) {
                
            m_cfg = cfg;                        
            m_dist.localFrame   = Frame ( m_dummyHit.normal );
            m_dist.storedLobes  = (int) m_cfg.fitting.nComponentsEnviroDirections;
            m_dist.randomInit( m_rnd );

            if ( m_particles.size() == 0 ) {
                rasterize( m_dist, enviro, raster, debug );
                return;
            }            
            fitting( m_particles, m_dist, /*true*/false );
            rasterize( m_dist, enviro, raster, debug );
        }


        IMPORTANCE_INLINE void refresh( const EnviroMap * enviro, 
            ImBitmap<Float> & raster,  DebuggingBundle & debug ) {            
            if ( m_particles.size() == 0 ) {
                rasterize( m_dist, enviro, raster, debug );
                return;
            }
            fitting( m_particles, m_dist, false );
            rasterize( m_dist, enviro, raster, debug );
        }


        IMPORTANCE_INLINE std::string toString() const {
            std::ostringstream ostr;
            ostr << "ParametricDistribution (enviro directions): [" << std::endl
                 << m_dist.toString() << std::endl
                 << "]" << std::endl;
            return ostr.str();
        }

    protected:
        IMPORTANCE_INLINE void fitting( const IStack<Particle>& particles, 
            GaussianMixtureSSESphere & dist, bool allowFirstFitStatic = false ) {

            OnlineEM::Info info;
            StepwiseEMConfig swcfg( m_cfg.fitting.alpha, m_cfg );            

            size_t count = particles.size();
            IIterator<const Particle> particlesIt = particles.cbegin();
            for ( size_t j = 0; j < count; j += BATCH_SIZE ) {
                Array samples( dist.localFrame, m_dummyHit, BATCH_SIZE );
                IIterator<const Particle> particlesEnd = 
                    ( particlesIt + BATCH_SIZE < particles.cend() ) ? 
                    particlesIt + BATCH_SIZE : particles.cend();
                fill( particlesIt, particlesEnd, samples );
                if ( j == 0 && allowFirstFitStatic ) {
                    OnlineEM::stepwiseStaticEM< GaussianMixtureModel, Array, QuadVector2, Float4, DIMENSION >
                        ( samples, (int) samples.size(), GaussianMixtureModel( dist ), 
                        QuadVector2(), m_stats, info, swcfg );
                } else {
                    OnlineEM::stepwiseEM< GaussianMixtureModel, Array, QuadVector2, Float4, DIMENSION >
                        ( samples, (int) samples.size(), GaussianMixtureModel( dist ), QuadVector2(), swcfg );
                }
            }
        }

        IMPORTANCE_INLINE void fill( IIterator<const Particle> & begin, IIterator<const Particle> & end, 
            Array & samples ) {            
            for ( auto it = begin; it < end; ++it ) {
                Particle p = *it;                
                p.incidentDir = p.position;
                p.weight = 1.f;
#ifdef LIBIMP_DEBUG
                p.position = Vector3( FLOAT_NAN );
#endif
                samples.add( p, m_dummyHit.normal, p.weight );
            }
            DataFilter<Array> filter( NULL, Config() );
            filter.clampAndStats( samples, &m_emInfo, m_cacheStats, m_cfg.fitting.enviroWeightClampFactor );
        }


        IMPORTANCE_INLINE void rasterize( const GaussianMixtureSSESphere & distribution, 
            const EnviroMap * enviro, ImBitmap<Float> & rastr, DebuggingBundle & debug ) const {
            const Vector2 dims( (Float) rastr.getWidth(), (Float) rastr.getHeight() );
            for ( int y = 0; y < rastr.getHeight(); ++y ) {
                for ( int x = 0; x < rastr.getWidth(); ++x ) {
                    const Vector3 dir   = EnviroSamplerMapping::getDir( Vector2( x + 0.5f, y + 0.5f ) / dims );
                    Float enviroPdf     = enviro != NULL ? enviro->_dirPdf( dir ) : 1.f;
                    Float impPdf        = distribution.pdf( dir );
                    rastr( x, y )       = enviroPdf * impPdf;
                    debug.m_pdf( x, y )  = impPdf;
                }
            }

            debug.updateHistogramAndCounts( m_particles );
        }


    private:
        void operator=( const ParametricDirections & ) {}
        GaussianMixtureSSESphere m_dist;
        Config m_cfg;        
        Hit m_dummyHit;
        const IStack<Particle>& m_particles;
        Random m_rnd;
        EMStats m_stats;

        EmInfo m_emInfo;
        CacheStats m_cacheStats;
    };
}