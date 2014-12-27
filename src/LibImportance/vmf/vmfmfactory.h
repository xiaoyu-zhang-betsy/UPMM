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


#if !defined( __VMFM_FACTORY_H )
#define __VMFM_FACTORY_H

#include <iostream>

#include "../shared/Config.h"
#include "../LibImportanceTypes.h"
#include "../shared/Particle.h"
#include "../shared/basicfactory.h"

#include "Vmf.h"
#include "datasetrecord.h"
#include "vmfminitializer.h"
#include "vmftypes.h"

#ifndef WIN32
#define WIN32
#include "Timer.h"
#undef WIN32
#else
#include "Timer.h"
#endif


#pragma warning(push)
#pragma warning(disable:4239)
#pragma warning(disable:4100)

namespace Importance {
    class VMFM;    


    class VMFMFactory : public BasicDistributionFactory {
    public:
        VMFMFactory() : 
            m_maxKappa( std::numeric_limits<Float>::max() ),
            m_timerEM( new Importance::Timer() ),
            m_timerTotal( new Importance::Timer() ) {}

        ~VMFMFactory() {
            delete m_timerEM;
			delete m_timerTotal;
        }

        // uloz si ty photony, budes je potrebovat pozdeji, protoze fit dava jen indexy do tohohle pole
        void init( const IStack<Importance::Particle>& photons, const Importance::Config * config );

        // Run E-M, store result to output and gather data about the E-M process
        virtual bool getInstance( 
            const Importance::Hit&      hit, 
            const KdQueryResult*        particles,
            int                         nParticlesFound, 
            DefaultDistributionModel&   output );

        virtual bool getInstanceProgressive( 
            const Importance::Hit   & hit,
            const KdQueryResult     * particles,
            int                     nParticlesFound,
            DefaultDistributionModel & output );

        bool runEM( 
            const VmfTypes::Dataset          &  dataset, 
            int                                 nDataset,
            Float                               wPhotonsUsed,
            size_t                              maxIter,
            VMFM                             &  output,
            Importance::EmInfo               &  info );


        virtual std::string toString() const {
            /// Check whether init was called first
            IMPORTANCE_ASSERT( m_photonMap != NULL );

            std::ostringstream oss;
            oss << "VMFMFactory[" << std::endl
                << "  nPhotonsPerComponent = " << m_config->vmfm.nPhotonsPerComponent << std::endl
                << "  maxEMIter = " <<  m_config->fitting.nMaxEMIter << std::endl                                      
                << "  m_nNNQuery = " << m_stats.nNNQuery << "," << std::endl
                << "  m_nPhotonsFoundStats = " << m_stats.nPhotonsFoundStats << "," << std::endl
                << "  avg. # photonsPerQuerry = " <<  
                    (m_stats.nNNQuery > 0 ? (m_stats.nPhotonsFoundStats / (double) m_stats.nNNQuery) : 0.0) 
                    << "," << std::endl
                << "  avg. radiusPerQuery = " << 
                    (m_stats.nNNQuery > 0 ? (m_stats.avgRadius / (double) m_stats.nNNQuery) : 0.0) << std::endl
                << "  avg. kapa = " << (m_stats.nNNQuery > 0 ? m_stats.avgKappa / (double)m_stats.nNNQuery : 0.0)
                    << std::endl
                << "  highest kapa = " << m_stats.highestKappa << std::endl
                << "  lowest kapa = " << m_stats.lowestKappa << std::endl
                << "  number of EM runs = " << m_stats.nEMRuns << std::endl
                << "  number of progressive EM runs" << m_stats.nProgressiveEMRuns << std::endl 
                << "  avg. number of EM iterations = " << 
                    (m_stats.nEMRuns > 0 ? m_stats.nEMIterations / (double) m_stats.nEMRuns : 0.0 ) << std::endl
				<< " total running time (ms) = " << 
                    (m_stats.runningTimeTotal / 1000) << std::endl
				 << " EM running time (ms) = " << 
                    (m_stats.runningTimeEM/ 1000 ) << std::endl
				<< "  avg. maximum likelihood " <<  (m_stats.nEMRuns > 0 ? m_stats.avgMaxLikelihood / (double) m_stats.nEMRuns : 0.0 ) << std::endl
                << "  # clamps = " << m_stats.nClamps << std::endl                 
                << "  # construction fail = " << m_stats.nConstructionFail << std::endl
                << "  # of trashed components = " << m_stats.nDisposedComponents << std::endl
                << "  # of re-normalizations = " << m_stats.nRenormalizations << std::endl
                << "]";
            return oss.str();
        }

    private:
        inline Float logLikelihood( const VMFM & dist, const VmfTypes::Dataset & dataset, int nDataset ) const;        

    private:
        struct ParamsRecord {
            ParamsRecord() : w( 0.f ), mi( 0.f ) {}

            Vector3 mi;
            Float w;
        };

        struct stats {
            Float avgRadius;
            Float highestKappa;
            Float lowestKappa;
            Float avgKappa;
            size_t nConstructionFail;
            size_t nEMRuns;
            size_t nClamps;            
            size_t nEMIterations;
            __int64 nPhotonsFoundStats;
            __int64 nNNQuery;
            size_t nDisposedComponents;
            size_t nRenormalizations;
            size_t nProgressiveEMRuns;
			
			Float avgMaxLikelihood;
			size_t runningTimeTotal;
			size_t runningTimeEM;
            struct stats() {
                avgRadius               = 0.0f;
                highestKappa            = std::numeric_limits<Float>::min();
                lowestKappa             = std::numeric_limits<Float>::max();
                avgKappa                = 0.0f;
                nConstructionFail       = 0;
                nEMRuns                 = 0;
                nClamps                 = 0;                
                nEMIterations           = 0;
                nPhotonsFoundStats      = 0;
                nNNQuery                = 0;
                nDisposedComponents     = 0;
                nRenormalizations       = 0;
                nProgressiveEMRuns      = 0;

				avgMaxLikelihood		= 0;
				runningTimeTotal		= 0;
				runningTimeEM			= 0;
            }
        } m_stats;

        VmfmInitializer m_vmfmInitializer;
        Float m_maxKappa;
        Float m_alpha;
    protected:
        
        Importance::Timer * m_timerEM;       
        Importance::Timer * m_timerTotal;
    };

}

#pragma warning(pop)

#endif