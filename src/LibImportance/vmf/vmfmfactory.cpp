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


#include <iostream>
#include <sstream>

#include "vmfmfactory.h"
#include "../shared/Utils.h"
#include "../shared/simplelogger.h"
#include "../em/stepwise_em.h"
#include "vmfm_stepwise_em.h"

#pragma warning(push)
#pragma warning(disable:4127)
#pragma warning(disable:4239)

namespace Importance {
    
    using namespace VmfTypes;

    void VMFMFactory::init( const IStack<Importance::Particle>& photons, const Importance::Config * config ) {
        BasicDistributionFactory::init( photons, config );
        m_maxKappa              = m_config->vmfm.maxKappa >= 0.f ? m_config->vmfm.maxKappa : VMFM::MAX_KAPA;
        m_alpha                 = clamp( m_config->fitting.alpha, 0.f, 1.f );
        
        if ( m_alpha > 0.f && m_alpha <= 0.5f ) {
            IMPORTANCE_ASSERT( false );
            ILog( EWarn, "Forgetting parameter alpha does not satisfy convergence conditions!" );
        }
    }


#if 0
    bool VMFMFactory::getInstance( 
        const Importance::Hit & hit, 
        const KdQueryResult* particles,
        const int nParticlesFound, 
        VMFM & output ) {

            m_stats.nNNQuery++;       
            IMPORTANCE_ASSERT( m_photonMap != NULL );
            IMPORTANCE_ASSERT( m_config->vmfm.nPhotonsPerComponent > 0 );

            /// probably wrong user settings
            if ( m_config->vmfm.nPhotonsPerComponent == 0 ) {
                m_stats.nConstructionFail++;                               
                output.m_info.constructionFail  = true;
                return false;
            }

            Dataset dataset;
            DatasetAdapter datasetAdapter( dataset );

            /* Pull out N-nearest particles from the map and store them into Dataset array  */
            Importance::DataFilter<DatasetAdapter> filter( m_photonMap, m_config->fitting.isFitOverSphere );
            if ( m_config->fitting.isUseStats ) {
                filter.setInfo( &output.m_info );
            }            
            filter.prepareData( particles, hit, nParticlesFound, datasetAdapter );
            int nParticlesUsed          = datasetAdapter.getNParticlesUsed();
            Float totalParticlesWeight  = datasetAdapter.getTotalWeight();

            /// there no data for fitting
            if ( nParticlesUsed == 0 ) {                                
                m_stats.nConstructionFail++;                               
                output.m_info.constructionFail  = true;
                return false;
            }        

            /// there are too few data for fitting (but return a distribution so that a cache record 
            /// could be set anyway)
            if ( nParticlesUsed < VMFM::MINIMAL_DATASIZE ) {
                output.initAsUninformed( hit.normal, 
                    harmonicDistanceFromInverseDistances( &dataset[0], &dataset[ nParticlesUsed ] ) );                
                output.init();                
                //IMPORTANCE_ASSERT(output.c[0].m_mi.isReal());
                return true;
            }
            
            /* Magicaly estimate the initial number of components based on the nn particles. */
            size_t nComponents = std::max( nParticlesUsed
                / m_config->vmfm.nPhotonsPerComponent, (size_t) 1 );
            nComponents = std::min( nComponents, m_config->vmfm.nMaxComponents );
            output = VMFM( nComponents );
            m_vmfmInitializer.init( m_config->vmfm.init, dataset, nParticlesUsed, hit.normal, output );
            output.init();
            output.k = 0;
            return true;
    }
#endif


    /// Returns true in case of success
    bool VMFMFactory::getInstance( 
                    const Importance::Hit & hit, 
                    const KdQueryResult* particles,
                    const int nParticlesFound, 
                    DefaultDistributionModel & output ) {

            IMPORTANCE_ASSERT( &dynamic_cast<VMFM&>(output) );
            VMFM & vmfmDist = (VMFM&) output;


            m_timerTotal->reset();
            m_stats.nNNQuery++;       
            IMPORTANCE_ASSERT( m_photonMap != NULL );
            IMPORTANCE_ASSERT( m_config->vmfm.nPhotonsPerComponent > 0 );

            /// probably wrong user settings
            if ( m_config->vmfm.nPhotonsPerComponent == 0 ) {
                m_stats.nConstructionFail++;                               
                vmfmDist.m_info.constructionFail  = true;
                return false;
            }

            Dataset dataset;
            DatasetAdapter datasetAdapter( dataset );
		                           
            /* Pull out N-nearest particles from the map and store them into Dataset array  */
            Importance::DataFilter<DatasetAdapter> filter( m_photonMap, m_config->fitting.isFitOverSphere );
            if ( m_config->fitting.isUseStats ) {
                filter.setInfo( &vmfmDist.m_info );
            }            
            filter.prepareData( particles, hit, nParticlesFound, datasetAdapter );
            int nParticlesUsed          = datasetAdapter.getNParticlesUsed();
            Float totalParticlesWeight  = datasetAdapter.getTotalWeight();

            /// there no data for fitting
            if ( nParticlesUsed == 0 ) {                                
                m_stats.nConstructionFail++;                               
                vmfmDist.m_info.constructionFail  = true;
                return false;
            }        

            /// there are too few data for fitting (but return a distribution so that a cache record 
            /// could be set anyway)
            if ( nParticlesUsed < VMFM::MINIMAL_DATASIZE ) {
                vmfmDist.initAsUninformed( hit.normal, 
                    harmonicDistanceFromInverseDistances( &dataset[0], &dataset[ nParticlesUsed ] ) );                
                vmfmDist.init();                
                //IMPORTANCE_ASSERT(vmfmDist.c[0].m_mi.isReal());
                return true;
            }

            /* Magicaly estimate the initial number of components based on the nn particles. */
            size_t nComponents = std::max( nParticlesUsed
                                            / m_config->vmfm.nPhotonsPerComponent, (size_t) 1 );
			nComponents = std::min( nComponents, m_config->vmfm.nMaxComponents );
            vmfmDist = VMFM( nComponents );

            bool success = false;
            const size_t N_RANDOM_RESTART = m_config->vmfm.nRandomRestarts;
            Float bestFitLkh = -std::numeric_limits<Float>::max();
            for ( size_t i = 0; i < N_RANDOM_RESTART; ++i ) {
                /* Keeps the currently supposed data model in particular EM iteration */
                VMFM curDist( nComponents );        

                //data model initalization for EM alg.  
                //m_timerInit->reset();
                m_vmfmInitializer.init( m_config->vmfm.init, dataset, nParticlesUsed, hit.normal, curDist );
                //vmfmDist.m_info.timeInit = m_timerInit->getMicroseconds();

#if defined(LIBIMP_DEBUG)
                Float error, weightsSum; bool isKappasValid;
                IMPORTANCE_ASSERT( curDist.check( error, isKappasValid, weightsSum ) );
#endif

				m_timerEM->reset();
                bool isDistValid = runEM( dataset, 
                                          nParticlesUsed, 
                                          totalParticlesWeight, 
                                          m_config->fitting.nMaxEMIter, 
                                          curDist, 
                                          vmfmDist.m_info );
				m_stats.runningTimeEM += m_timerEM->getMicroseconds();
                size_t last = std::max( (size_t) 0, vmfmDist.m_info.lkhFunction.size() - 1 );
                if ( isDistValid && last < vmfmDist.m_info.lkhFunction.size() && vmfmDist.m_info.lkhFunction[ last ] > bestFitLkh ) {
                    success = true;
                    vmfmDist = curDist;
                    bestFitLkh = vmfmDist.m_info.lkhFunction[ last ];
                }
            }

            if ( vmfmDist.nComponents() == 0 ) {
                /// if fitting happens to throw away all components we do not wan't to have an empty distribution
                vmfmDist.initAsUninformed( hit.normal, std::sqrt(1.f/dataset[0].getInverseDistance()) );
                ILog( EWarn, "All components were disposed during fitting. \
                             Initializing distribution as uninformed." );
                success = true;
            }

            if ( !success ) {
                return false;
            }

            /* Prepare for sampling */            
            vmfmDist.init();

            /// TODO: Importance Lib
		    //m_mut->lock();
			    m_stats.avgKappa += vmfmDist.getAvgKappa();
			    m_stats.nPhotonsFoundStats += nParticlesUsed;			    
                m_stats.highestKappa = std::max( m_stats.highestKappa, vmfmDist.getMaxKappa() );
                m_stats.lowestKappa = std::min( m_stats.lowestKappa, vmfmDist.getMinKappa() );
		    //m_mut->unlock();
            m_stats.runningTimeTotal += m_timerTotal->getMicroseconds();
            vmfmDist.k = (unsigned int) nParticlesUsed; //remember first batch (0 means forget statistics for stepwiseEM)            
            return nParticlesUsed > 0;    
    }

    bool VMFMFactory::getInstanceProgressive( 
        const Importance::Hit & hit,
        const KdQueryResult* particles,
        const int nParticlesFound,
        DefaultDistributionModel& output) {

        IMPORTANCE_ASSERT( &dynamic_cast<VMFM&>(output) );
        VMFM & vmfmDist = (VMFM&) output;

        m_stats.nProgressiveEMRuns++;
        IMPORTANCE_ASSERT( vmfmDist.nComponents() <= Importance::MAX_FITTED_LOBES );
        IMPORTANCE_ASSERT( vmfmDist.nComponents() > 0 );
        if ( nParticlesFound < 1 || vmfmDist.nComponents() == 0 ) {
            std::stringstream msg;
            msg << "Zero photons or zero # of components. nPhotonsFound = " <<  nParticlesFound
                << "; nComponents = vmfmDist.nComponents()";
            ILog( EWarn, msg.str() );
            return false;
        }   
        IMPORTANCE_ASSERT( m_photonMap != NULL );

        const size_t nComponents = vmfmDist.nComponents();


        Dataset dataset;
        DatasetAdapter datasetAdapter( dataset );

        /* Pull out N-nearest particles from the map and store them into Dataset array  */
        Importance::DataFilter<DatasetAdapter> filter( m_photonMap, m_config->fitting.isFitOverSphere );
        if ( m_config->fitting.isUseStats ) {
            filter.setInfo( &vmfmDist.m_info );
        }            
        filter.prepareData( particles, hit, nParticlesFound, datasetAdapter );
        int nParticlesUsed          = datasetAdapter.getNParticlesUsed();
        //Float totalParticlesWeight  = datasetAdapter.getTotalWeight();


        if ( nParticlesUsed == 0 ) {
            return false;
        }
    
        StepwiseEMConfig swcfg( m_alpha, *m_config );
        OnlineEM::stepwiseEM( dataset, nParticlesUsed, EM::VmfmMixtureModel( vmfmDist, m_maxKappa ), Vector3(), swcfg );

        Float error         = 0.0f;
        bool isKappasValid  = true;
        Float weightsSum    = 0.f;
        bool isValid        = vmfmDist.check( error, isKappasValid, weightsSum );

        if ( !isValid && isKappasValid && weightsSum > 0.f ) {
            /// re-normalize
            for ( size_t h = 0; h < nComponents; ++h ) {
                vmfmDist.w[ h ] /= weightsSum;
            }
            m_stats.nRenormalizations++;            
            ILog( EWarn, "Distribution check failed - weights" );
            return false;
        } else if ( !isValid ) {
            /// Some of kappas is/are in invalid range or weights sum to zero
            IMPORTANCE_ASSERT( false );
            ILog( EWarn, "Distribution check failed - kappas" );
            m_stats.nConstructionFail++;            
            return false;
        }
        
        vmfmDist.init();

        return true;
    }

    bool VMFMFactory::runEM( const Dataset & dataset, 
                             int nDataset, 
                             Float wPhotonsUsed,
                             size_t maxIter,                              
                             VMFM & output,
                             Importance::EmInfo & /*info*/ ) {
        const size_t nComponents    = output.nComponents();
        size_t iIter                = 0;
        bool stop                   = false;
        size_t lastClampAt          = 0;
        bool warningFlag            = false;
        Float maxLkhEst             = -std::numeric_limits<Float>::max();
        Float lkhFunctionMax        = -std::numeric_limits<Float>::max();
        Float prevMaxLkhEst         = -std::numeric_limits<Float>::max();
        
        /* Posterior pdf of sampling i-th point from h-th component */
        IStack<Float> pst( nComponents, 0.0f );
        /* Temporary storage for computed values */
        VMFM newDist( nComponents );
        VMFM dist( output );
        /* Temporary memory used for computation of harmonic distances */
        IStack<Float> tmpHarmDist( nComponents, 0.f );
        
        //for debugging only
        IStack<Float> lkhFunction;                           
        size_t nPhotonsIgnoredTotal = 0;


        for ( ; iIter < maxIter && !stop; ++iIter ) {  
            //std::cout << "\r" << iIter << "/" << maxIter;

            prevMaxLkhEst           = maxLkhEst;
            maxLkhEst               = 0.0f;
            Float wPhotonsIgnored   = 0.0f;
			size_t nPhotonsIgnored	= 0;

            std::fill(tmpHarmDist.begin(), tmpHarmDist.begin()+nComponents, 0.f);

            //iterate through all points            
            for ( int i = 0; i < nDataset; ++i ) {
                Float denom = 0.0f;
                /* Weight of direction set according to corresponding photon value */
                Float dirWeight = dataset[ i ].value;
                
                //compute posterior pdf for point i and normalization term 
                for ( size_t h = 0; h < nComponents; ++h ) {
                    pst[ h ] = dist.w[ h ] * dist.c[ h ].pdf( dataset[ i ].dir );
                    denom += pst[ h ];
                }

                if ( denom < std::numeric_limits<Float>::min() /*1e-36f*/ ) {
                    /// we ignore the photon as it is too "far" from every lobe                    
                    wPhotonsIgnored += dirWeight;
					nPhotonsIgnored++;
                    continue;
                }

                //compute weights and mis from posterior probabilities
                for ( size_t h = 0; h < nComponents; ++h ) {
                    /// Do not pre-compute the division for numerical stability reasons!
                    Float posterior = dirWeight * pst[ h ] / denom;
                    
                    IMPORTANCE_ASSERT( !( Importance::isNaN( posterior ) ) && 
                                       !( Importance::isINF( posterior ) ) );

                    newDist.w[ h ] += posterior;
                    newDist.c[ h ].m_mi += dataset[ i ].dir * posterior;
                        
                    /// use this later for normalization
                    newDist.c[ h ].m_harmonicDistance += posterior;
                    /// harmonic distance sum
                    tmpHarmDist[ h ] += posterior * dataset[ i ].invDist;
                    //IMPORTANCE_ASSERT( i < nDataset - 1 || (!Importance::isNaN( tmpHarmDist[ h ] ) && !Importance::isINF( tmpHarmDist[ h ] )) );
                    //IMPORTANCE_ASSERT( i < nDataset - 1 || (tmpHarmDist[ h ] >= std::numeric_limits<Float>::min()) );
                }

                maxLkhEst += dirWeight * std::log( denom );
            }

			nPhotonsIgnoredTotal += nPhotonsIgnored;

            lkhFunction.push( maxLkhEst );

            /// stopping condition - if the value drops twice in a row under so far reached maximum of lkhFunction
            /// or the improvement was twice not high enough
            //bool isUnderThr = maxLkhEst - prevMaxLkhEst < VMFM::MINIMAL_STEP;
            //stop            = warningFlag && ( lkhFunctionMax >= maxLkhEst || isUnderThr );                                   
            //warningFlag     = ( lkhFunctionMax >= maxLkhEst || isUnderThr );

            bool isUnderThr = std::abs( prevMaxLkhEst / maxLkhEst - 1.f ) < VMFM::MINIMAL_ERROR;
            stop            = warningFlag && isUnderThr;                                   
            warningFlag     = isUnderThr;

            lkhFunctionMax  = std::max( lkhFunctionMax, maxLkhEst );

            if ( stop ) {                                
                break;
            }            

            // debuging - leave this peace of code later
            dist.invalidate();
            //iterate per component, finish M-step, clear newDist for next iteration
            for ( size_t h = 0; h < nComponents; ++h ) {
                // ignore components which has already zero weights
                if ( newDist.w[ h ] <= 1e-7 ) {
                    dist.w[ h ] = 0.f;
                    //dist.c[ h ].m_kappa = 0.f;
                    dist.c[ h ].initPDFOnly();
                    continue;
                }

                Float length            = newDist.c[ h ].m_mi.length();
                length                  = std::max( length, std::numeric_limits<Float>::min() );                
                /// When length comes close to zero (hence the max function) it means, that no photon is likely to come from 
                /// the h-th lobe (SHOULDN'T WE DISCARD THE COMPONENT?)
                /// This might happened with caustics
                dist.c[ h ].statsMi  = newDist.c[ h ].m_mi;
                dist.c[ h ].m_mi     = newDist.c[ h ].m_mi / length;
                Float r              = length / newDist.w[ h ];
                IMPORTANCE_ASSERT( r > 0.f );
                    
                /// MAGIC clamping Value - can be derived from formula for kappa (we mustn't devide by zero)  
                /// TODO: derive it for double precision as well!
                //r                       = std::min( r, 0.999994f );
#ifdef IMPORTANCE_SINGLE_PRECISION
                r                       = std::min( r, Float(0.99999997f) );
#else
                /// TODO: check this for double - hasn't been used yet
                r                       = std::min( r, Float(0.999999999997f) );
#endif
                Float rSqr              = r * r;
                dist.c[ h ].m_kappa     = r * ( 3.0f - rSqr ) / ( 1.0f - rSqr );
                dist.w[ h ]             = newDist.w[ h ] / ( wPhotonsUsed - wPhotonsIgnored );
                dist.c[ h ].m_harmonicDistance = newDist.c[ h ].m_harmonicDistance;

                /// for online version
                dist.c[ h ].statsMi     /=  ( wPhotonsUsed - wPhotonsIgnored );
                dist.c[ h ].statsW      = dist.w[ h ];         

                /// !!If this happened regularly in some scenes then clamping by MAX_KAPA is not 
                /// enough and you should solve this somehow better!!
                /// This is severe problem as you have just lost the most concentrated lobe
                //if ( Importance::isNaN(dist.c[ h ].m_kappa) || Importance::isINF(dist.c[ h ].m_kappa) || dist.c[ h ].m_kappa < 0.f ) {
                //    std::cout << "(runEM)Input dist: " << output.toString() << std::endl
                //        << "kappa: " << dist.c[ h ].m_kappa << "; r = " << r << "; dist.c[ h ].wRaw = " << dist.c[ h ].wRaw
                //        << std::endl << "length = " << length << ";newDist.w[ h ] = " << newDist.w[ h ] << std::endl; 
                //}
                
                // TODO: SOLVE THIS BETTER! INVESTIGATE WHY IS IT HAPPENING FROM TIME TO TIME IN THE FIRST PLACE
                IMPORTANCE_ASSERT( dist.c[ h ].m_kappa >= 0.f );
                IMPORTANCE_ASSERT( !Importance::isNaN( dist.c[ h ].m_kappa ) && !Importance::isINF( dist.c[ h ].m_kappa ) );  

                /// Clamp kappa
                if ( dist.c[ h ].m_kappa < 0.f || 
                        dist.c[ h ].m_kappa > m_maxKappa ) {
                    /// Continue in fitting but keep the kapa under MAX_KAPA so that we 
                    /// won't get into the trouble                        
                    dist.c[ h ].m_kappa = clamp( dist.c[ h ].m_kappa, 0.f, m_maxKappa );
                    lastClampAt = iIter;
                    /// TODO: Importance Lib
                    //m_mut->lock();
                        ++m_stats.nClamps;
                    //m_mut->unlock();
                }
                
                //precompute constants (used for sampling and pdf computation) for h-th lobe
                dist.c[ h ].initPDFOnly();


                //clear newDist for next iteration
                newDist.w[ h ] = 0.0f;
                newDist.c[ h ].m_mi = Vector3( 0.0f );
                newDist.c[ h ].m_harmonicDistance = 0.f;                    
            }
            
            Float error         = 0.0f;
            bool isKappasValid  = true;
            Float weightsSum    = 0.f;
            bool isValid        = dist.check( error, isKappasValid, weightsSum );
            
            if ( !isValid && isKappasValid && weightsSum > 0.f ) {
                /// re-normalize
                for ( size_t h = 0; h < nComponents; ++h ) {
                    dist.w[ h ] /= weightsSum;
                }
                m_stats.nRenormalizations++;
                ILog( EWarn, "Distribution check failed." );
            } else if ( !isValid ) {
                /// Some of kappas is/are in invalid range or weights sum to zero
                IMPORTANCE_ASSERT( false );
                m_stats.nConstructionFail++;
                output.m_info.constructionFail  = true;
                return false;
            }

        } // end of EM iteration

        /// if EM stopped because of maximum number of iterations limit, compute log-likelihood of current fit
        if ( !stop ) {
            lkhFunction.push( logLikelihood( dist, dataset, nDataset ) );
        }

        /* debug staff */
        output.m_info.lkhFunction        = lkhFunction;
        output.m_info.nIteration         = iIter;
        /// carefull, one photon could have been ignored multiple times during EM
        output.m_info.nPhotonsIgnored    = nPhotonsIgnoredTotal;

        /// Copy only components with non-zero weight
        const Float minWeight = std::numeric_limits<Float>::min();                    
        output = VMFM();
        output.reserve( nComponents );
        size_t k = 0;
        for ( size_t h = 0; h < nComponents; ++h ) {
            if ( dist.w[ h ] > minWeight ) {
                output.c.push( dist.c[ h ] );
                output.w.push( dist.w[ h ] );                    
                /* Compute harmonic distance - used for interpolation */
                output.c[ k ].m_harmonicDistance = dist.c[ h ].m_harmonicDistance / tmpHarmDist[ h ];   
                /// TODO: really? solve this better
                if ( output.c[ k ].m_harmonicDistance > std::numeric_limits<Float>::max() ) {
                    output.c[ k ].m_harmonicDistance = 1.f / dataset[ 0 ].invDist;
                }
                //if ( Importance::isNaN( output.c[ k ].m_harmonicDistance ) || 
                //     Importance::isINF( output.c[ k ].m_harmonicDistance ) ) {
                //Log ( EInfo, "PROBLEM: %s, tmpHarmDist[ h ] = %e \n output: %s, k = %u", dist.toString().c_str(), tmpHarmDist[ h ], output.toString().c_str(), k );
                //}
                //IMPORTANCE_ASSERT( !Importance::isNaN( output.c[ k ].m_harmonicDistance ) && 
                //             !Importance::isINF( output.c[ k ].m_harmonicDistance ) );
                k++;
            } else {
                std::stringstream msg;
                msg << "Component " << h << " was disposed: " << dist.c[ h ].toString().c_str() << " , w[ " << h << " ] = " << dist.w[ h ];
                ILog( EInfo, msg.str() );
                m_stats.nDisposedComponents++;
            }
        }

        if ( k == 0 ) {            
            ILog( EWarn, "Construction failed: all components had zero weight." );
            return false;
        }
        output.m_nComponents = k;

        m_stats.nEMRuns++;        
        m_stats.nEMIterations += iIter;
		m_stats.avgMaxLikelihood += maxLkhEst;
        return true;
    }


    Float VMFMFactory::logLikelihood( const VMFM & dist, const Dataset & dataset, int nDataset ) const {
        Float res = 0.f;
        for ( int i = 0; i < nDataset; ++i ) {
            Float tmp = 0.0f;
            Float dirWeight = dataset[ i ].value;
            //compute posterior pdf for point i and normalization term 
            for ( size_t h = 0; h < dist.nComponents(); ++h ) {
                tmp += dist.w[ h ] * dist.c[ h ].pdf( dataset[ i ].dir );                
            }

            res += dirWeight * std::log( tmp );
        }

        return res;
    }

}

#pragma warning(pop)