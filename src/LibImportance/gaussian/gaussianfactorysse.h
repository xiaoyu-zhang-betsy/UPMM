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

#pragma warning(push)
#pragma warning(disable:4127)
#pragma warning(disable:4239)

#include "../shared/Config.h"
#include "../shared/basicfactory.h"
#ifndef WIN32
    #define WIN32
#endif
#include "../vmf/Timer.h"
#include "../em/stats.h"

#include "arraysse.h"
#include "gaussiansse.h"
#include "../em/stepwise_emsse.h"
#include "gaussian_stepwise_emsse.h"

namespace Importance {    
	template<typename TScalar, int DIMENSION, typename TDistribution>
    class GaussianMixtureFactorySIMD : public BasicDistributionFactory {
    public:
//		typedef typename GuassianMultiLobe< TScalar >::QuadVector2 QuadVector2;
		typedef EM::GaussianMixtureModel< TScalar, DIMENSION, TDistribution > GaussianMixtureModel;        
		typedef ArrayGaussianSSE< TScalar, DIMENSION, TDistribution > Array;

		GaussianMixtureFactorySIMD() :
			m_timerEM( new Importance::Timer() ),
			runningTimeEM( 0 ) {}

		~GaussianMixtureFactorySIMD() { delete m_timerEM; }

    /************************************************************************/
    /* BasicDistributionFactory iface */
    /************************************************************************/

        inline void init( const IStack<Importance::Particle>& particles, const Importance::Config * config ) {
            BasicDistributionFactory::init( particles, config );
            /// never try fitting to all data over sphere (we use mapping only for hemisphere)
            IMPORTANCE_ASSERT( m_config->fitting.isFitOverSphere == false );
            m_alpha = clamp( m_config->fitting.alpha, 0.f, 1.f );
            if ( m_alpha > 0.f && m_alpha <= 0.5f ) {
                IMPORTANCE_ASSERT( false );
                ILog( EWarn, "Forgetting parameter alpha does not satisfy convergence conditions!" );
            }
        }

	/************************************************************************/
    /* AbstractDistributionFactory interface */
    /************************************************************************/

        /** 
         * Use this method when new cache record 
         * Run E-M, store result to output and gather data about the E-M process 
         */
        virtual bool getInstance( 
            const Importance::Hit& hit, 
            const KdQueryResult* particles,
            int nParticlesFound, 
            DefaultDistributionModel & outDist) {
            
            TDistribution & gmm = gmm_cast( outDist );

            if ( m_config->fitting.isPureOnline ) {
                IMPORTANCE_ASSERT( gmm.nComponents() == 0 );

                gmm = TDistribution((int)m_config->fitting.nMaxComponents);
                gmm.localFrame = Frame(hit.normal);
                gmm.randomInit( m_rnd );
                return getInstanceProgressive( hit, particles, nParticlesFound, gmm );                
            }
                        
            return getInstanceStatic( hit, particles, nParticlesFound, gmm );
        }

		
		/** 
        * Use this method with existing cache records to their progressive update 
        * Run E-M, store result to output and gather data about the E-M process 
        */
		virtual bool getInstanceProgressive( 
			        const Importance::Hit& hit, 
			        const KdQueryResult* particles,
			        int nParticlesFound,                    
			        DefaultDistributionModel & outDist ) {                  
            TDistribution & gmm = gmm_cast( outDist );
            IMPORTANCE_ASSERT( gmm.nComponents() > 0 );
            IMPORTANCE_ASSERT( m_photonMap != NULL );
            IMPORTANCE_ASSERT( gmm.nComponents() <= Importance::MAX_FITTED_LOBES );

            const Frame & frame = gmm.localFrame;
            /* Pull out N-nearest particles from the map on the proper hemisphere and push them into samples array */            
            Importance::DataFilter<Array> filter( m_photonMap, *m_config );
            IMPORTANCE_INC_STATS( m_stats.nNNQuery );
            IMPORTANCE_INC_STATS_MORE( m_stats.nPhotonsFound, nParticlesFound );
            IMPORTANCE_INC_STATS_MORE( gmm.info.nPhotonsFound, nParticlesFound );
            Array samples( frame, hit,  nParticlesFound );
            filter.prepareData( particles, hit, nParticlesFound, gmm.m_cacheStats, samples );            
            filter.clampAndStats( samples, EMINFO( gmm ), gmm.m_cacheStats, m_config->fitting.weightClampFactor);
            
            /// there are no data for fitting
            if ( samples.size() == 0 ) {                                
                return false;
            }            
                   
            m_timerEM->reset();

            StepwiseEMConfig swcfg( m_alpha, *m_config );                
            OnlineEM::stepwiseEM< GaussianMixtureModel, Array, QuadVector2, TScalar, DIMENSION >
				( samples, (int) samples.size(), GaussianMixtureModel( gmm ), QuadVector2(), swcfg );

            gmm.setUndersampled( gmm.isUndersampled() && 
                gmm.m_cacheStats.photonCount < m_config->cache.undersampledCoefOnline * m_config->fitting.nMaxComponents );


            IMPORTANCE_INC_STATS_MORE( m_stats.runningTimeEM, m_timerEM->getMicroseconds() );
#ifdef LIBIMP_STATS
            IMPORTANCE_INC_STATS_MORE( gmm.info.nPhotonsUsed, samples.size() );
			if (m_config->fitting.isUseStats) {
                gmm.info.particles.reserve( samples.size() );
                for ( Array::iterator it = samples.begin(); it != samples.end(); ++it ) 
                {
                    Particle p;
                    p.incidentDir = -it->original;
                    p.weight = it->weight[ 0 ];
                    gmm.info.particles.push( p );
                }
            }
#endif            
            return true;
		}

    /************************************************************************/
    /* Print */
    /************************************************************************/

        virtual std::string toString() const {
            /// Check whether init was called first
            IMPORTANCE_ASSERT( m_photonMap != NULL );
            std::ostringstream oss;
            oss << "GaussianFactory[" << std::endl
                << "  nPhotonsPerComponent = " << m_config->vmfm.nPhotonsPerComponent << std::endl
				<< "  avg. components = " <<  
                    (m_stats.nNNQuery > 0 ? (m_stats.avgComponents / (double) m_stats.nNNQuery) : 0.0) << std::endl
                << "  maxEMIter = " <<  m_config->fitting.nMaxEMIter << std::endl                                      
                << "  m_nNNQuery = " << m_stats.nNNQuery << "," << std::endl
                << "  m_nPhotonsFoundStats = " << m_stats.nPhotonsFound << "," << std::endl
                << "  avg. # photonsPerQuerry = " <<  
                    (m_stats.nNNQuery > 0 ? (m_stats.nPhotonsFound / (double) m_stats.nNNQuery) : 0.0) 
                    << "," << std::endl
                << "  number of EM runs = " << m_stats.nEMRuns << std::endl
                << "  avg. number of EM iterations = " << 
                    (m_stats.nEMRuns > 0 ? m_stats.nEMIterations / (double) m_stats.nEMRuns : 0.0 ) << std::endl
				 << " EM running time (ms) = " << 
                    (m_stats.runningTimeEM/ 1000 ) << std::endl
				<< "  avg. maximum likelihood " <<  (m_stats.nEMRuns > 0 ? m_stats.avgMaxLikelihood / (double) m_stats.nEMRuns : 0.0 ) << std::endl
                << "  # construction fail = " << m_stats.nConstructionFail << std::endl
                << "  # of trashed components = " << m_stats.nDisposedComponents << std::endl
                << "]";
            return oss.str();
        }

		inline const EMStats & getStatistics() const {
			return m_stats;
		}
	protected:
		       
		bool getInstanceStatic( 
			const Importance::Hit& hit, 
			const KdQueryResult* particles,
			int nParticlesFound, 
			TDistribution & outDist ) {       
//            IMPORTANCE_ASSERT( outDist.nComponents() == 0 );
            IMPORTANCE_ASSERT( m_photonMap != NULL );                   
            
			Frame frame(hit.normal);           
            /* Pull out N-nearest particles from the map on the proper hemisphere and push them into samples array */            
            Importance::DataFilter<Array> filter( m_photonMap, *m_config );
			IMPORTANCE_INC_STATS( m_stats.nNNQuery );
            IMPORTANCE_INC_STATS_MORE( m_stats.nPhotonsFound, nParticlesFound );
            IMPORTANCE_INC_STATS_MORE( outDist.info.nPhotonsFound, nParticlesFound );
			Array samples( frame, hit, nParticlesFound );            
            filter.prepareData( particles, hit, nParticlesFound, outDist.m_cacheStats, samples );
            filter.clampAndStats( samples, EMINFO( outDist ), outDist.m_cacheStats, m_config->fitting.weightClampFactor);           
            
            size_t nComponents = m_config->fitting.nMaxComponents;
            if ( samples.size() == 0 || 
                 nComponents * m_config->vmfm.nPhotonsPerComponent > samples.size() ) {                
                outDist.localFrame = frame;
                outDist.storedLobes = (int)nComponents;
                /* Initialize the distribution based on the hit normal. */			
                outDist.randomInit( m_rnd );
                outDist.setUndersampled( true );                
                return true;
            }
            
			Float likelihood;
			m_timerEM->reset();

            StepwiseEMConfig swcfg( m_alpha, *m_config );                             
            OnlineEM::Info info;

            /* Random restarts */            
            likelihood = -std::numeric_limits<Float>::max();
            for ( int i = 0; i < (int) m_config->fitting.nRandomRestarts + 1; ++i ) { 
                TDistribution tmpDist = outDist;                
                tmpDist.localFrame = frame;
                tmpDist.storedLobes = (int)nComponents;
                tmpDist.randomInit( m_rnd );                
				OnlineEM::stepwiseStaticEM< GaussianMixtureModel, Array, QuadVector2, TScalar, DIMENSION >
					( samples, (int) samples.size(), GaussianMixtureModel( tmpDist ), QuadVector2(), m_stats, info, swcfg );
                if ( info.lkh > likelihood ) {
                    likelihood  = info.lkh;
                    outDist     = tmpDist;
                }
            }          
            outDist.setUndersampled( info.nAdjustedSamples < m_config->cache.undersampledCoef * m_config->fitting.nMaxComponents );

			IMPORTANCE_INC_STATS_MORE( m_stats.runningTimeEM, m_timerEM->getMicroseconds() );
            IMPORTANCE_INC_STATS_MORE( m_stats.avgMaxLikelihood, likelihood );
#ifdef LIBIMP_STATS
            IMPORTANCE_INC_STATS_MORE( outDist.info.nPhotonsUsed, samples.size() );
			if (m_config->fitting.isUseStats) {
                outDist.info.particles.reserve( samples.size() );
				for ( Array::iterator it = samples.begin(); it != samples.end(); ++it ) 
				{
					Particle p;
					p.incidentDir = -it->original;
					p.weight = it->weight[ 0 ];
					outDist.info.particles.push( p );
				}
			}
#endif
            return true;
		}

		/// piggy stuff but we have no time to polish the architecture now
        inline IMPORTANCE_INLINE TDistribution & gmm_cast( DefaultDistributionModel & dist ) const {
            IMPORTANCE_ASSERT( &dynamic_cast<TDistribution &>( dist ) );
            return (TDistribution &) dist;
        }

	private:
		Importance::Timer * m_timerEM;
		size_t runningTimeEM;
        Float m_alpha;
        Random m_rnd;
		EMStats m_stats;
    };
}

#pragma warning(pop)