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
#include "stats.h"
#include "stepwise_EM_config.h"

#undef ADJUST_DATASET
//#define ADJUST_DATASET 1
#define ADJUST_DATASET 0

namespace Importance {  

    namespace OnlineEM {
        const int MAX_DATASETSIZE = 4 * MAX_KNN_PARTICLES;
        const Float MIN_WEIGHT_OEM    = 0.f;

        template<typename TVector, typename TScalar>
        struct Sample {
            IMPORTANCE_INLINE const TVector & getPoint() const { return m_point; }
            IMPORTANCE_INLINE const TScalar & getValue() const { return m_weight; }
            IMPORTANCE_INLINE void setValue( const TScalar & w ) { m_weight = w; }

            TVector m_point;
            TScalar m_weight;
        };

        struct Info {
            Float lkh;
            int nIter;
            int nAdjustedSamples;
        };

        template<typename TDataSet, typename TVector, typename TScalar>
        IMPORTANCE_INLINE int adjustDataset( const TDataSet & dataset, int n, ImStaticArray<OnlineEM::Sample<TVector, TScalar>,MAX_DATASETSIZE> & samples, float minw, bool & warn ) {
            IMPORTANCE_ASSERT( n > 0 );
            
            warn             = false;

            // compute average value
            Float maxw       = 0.f;
            for ( int i = 0; i < n; ++ i ) {
                maxw += dataset[ i ].getValue()[ 0 ];
            }
            maxw /= n;

            int j = 0;
            for ( int i = 0; i < n; ++ i ) {
                Float q = dataset[ i ].getValue()[ 0 ];
                //IMPORTANCE_ASSERT( q > minw );
                if ( q > minw ) {                                                                
                    while ( q > maxw ) {
                        IMPORTANCE_ASSERT( j < MAX_DATASETSIZE );
                        q -= maxw;
                        samples[ j ].m_point           = dataset[ i ].getPoint();
                        samples[ j ].setValue( TScalar( maxw ) );
                        j ++;                        
                    }
                    if ( q > maxw * 1e-3f ) {
                        IMPORTANCE_ASSERT( j < MAX_DATASETSIZE );
                        samples[ j ].m_point           = dataset[ i ].getPoint();
                        samples[ j ].setValue( TScalar( q ) );
                        j ++;
                    }
                } else {
                    warn = true;
                }
            }

            std::random_shuffle( samples.begin(), samples.begin() + j );

            return j;
        }

        template<typename TMixtureModel, typename TDataSet, typename TVector, typename TScalar, int DIMENSION>
        void stepwiseEM( 
            const TDataSet & inputDataset,
            int nDatasetSize,                        
            TMixtureModel & model, 
            const TVector &,
            const StepwiseEMConfig & cfg = StepwiseEMConfig() ) { 
                IMPORTANCE_ASSERT( nDatasetSize > 0 );
                IMPORTANCE_ASSERT( cfg.alpha > 0.5f && cfg.alpha <= 1.f );            
           
#if ADJUST_DATASET                
                bool warn        = false;
                ImStaticArray<Sample<TVector, TScalar>,MAX_DATASETSIZE> dataset;
                nDatasetSize = adjustDataset<TDataSet,TVector,TScalar>( inputDataset, nDatasetSize, dataset, MIN_WEIGHT_OEM, warn );
                if ( warn ) {
                    ILog( EWarn, "There are weights smaller than %f!", MIN_WEIGHT_OEM );
                }
#else
                const TDataSet & dataset = inputDataset;
#endif

                ImStaticArray<TScalar,MAX_FITTED_LOBES / DIMENSION> pdf;
//                int batchSize = clamp( cfg.batchSize, 1, nDatasetSize );
				const int lobesCount = model.nComponents() / DIMENSION;
				static const TScalar v1 = TScalar( 1.0f );
                for ( int i = 0; i < nDatasetSize; ++i ) {

					//////////////////////////////////////////////////////////////////////////
					/// E-step: update sufficient statistics
					//////////////////////////////////////////////////////////////////////////

                    /// compute posterior probabilities                
                    std::fill( pdf.begin(), pdf.begin() + lobesCount, TScalar(0.f) );
                    TScalar denom				  = TScalar( 0.f );
                    const TScalar & dataWeight    = dataset[ i ].getValue();            
                    const TVector & x			  = dataset[ i ].getPoint();

                    for ( int h = 0; h < lobesCount; ++h ) {
                        pdf[ h ] = model.responsibility( h, x );
                        denom += pdf[ h ];
                    }
					denom = TScalar( TScalar::dot( denom, v1 ) );
                    auto stats  = model.getStats();
                    Float eta   = std::pow( stats.getCount() + 1.f, -cfg.alpha );
                    /// update weight statistic (common for all components)
                    stats.updateWeight( eta, dataWeight[ 0 ] );
					TScalar etaV( eta );

                    for ( int h = 0; h < lobesCount; ++h ) {
                        TScalar posterior = pdf[ h ] / denom;
                        stats.updateMultiLobe( h, etaV, x, dataWeight, posterior );
                    }            
                    stats.increment();

					//////////////////////////////////////////////////////////////////////////
					/// M-step: update model
					//////////////////////////////////////////////////////////////////////////         
                
					//if ( model.getStats().getCount() + 1 > cfg.delayedUpdate ) {
					if ( model.getStats().getCount() % cfg.batchSize == 0  ) {
						for ( int h = 0; h < lobesCount; ++h ) {
							model.update( h, cfg.regularizationType );
						}
					}

					IMPORTANCE_ASSERT( model.check() );
                    
                }								
        }



        /// General way to compute log-likelihood of a given mixture (Gaussian, vMFm, ...)
        template<typename TMixtureModel, typename TDataSet, typename TScalar, int DIMENSION >
        class MixtureLikelihoodFunctional {
            void operator=(const MixtureLikelihoodFunctional&);
        public:
            MixtureLikelihoodFunctional( const TMixtureModel & model ) 
                : m_model( model ) {}

            /// Compute log-likelihood for given dataset
            Float eval( const TDataSet & dataset, int nDatasetSize ) const {
				const int lobesCount = m_model.nComponents() / DIMENSION;
                Float res = 0;
                for ( int i = 0; i < nDatasetSize; ++i ) {

                    TScalar val = TScalar( 0.f );
                    for ( int h = 0; h < lobesCount; ++h ) {                 
                        val += m_model.responsibility( h, dataset[ i ].getPoint() );
                    }

                    res += dataset[ i ].getValue()[ 0 ] * ( Ff::log( TScalar::dot( val , TScalar( 1.0f ) ) ) );
                }

                return res;
            }

        private:
            const TMixtureModel & m_model;        
        };

        template<typename TMixtureModel, typename TDataSet, typename TVector, typename TScalar, int DIMENSION>
        void stepwiseStaticEM( 
            const TDataSet & inputDataset,
            int nDatasetSize,                        
            TMixtureModel & model, 
            const TVector &,
            EMStats & emStats,
            Info & info,       
            const StepwiseEMConfig & cfg = StepwiseEMConfig() ) { 
                IMPORTANCE_ASSERT( nDatasetSize > 0 );
                IMPORTANCE_ASSERT( cfg.alpha > 0.5f && cfg.alpha <= 1.f );         

#if ADJUST_DATASET
                typedef ImStaticArray<Sample<TVector, TScalar>,MAX_DATASETSIZE> DatasetType;
                DatasetType dataset;
                bool warn    = false;
                nDatasetSize = adjustDataset<TDataSet,TVector>( inputDataset, nDatasetSize, dataset, MIN_WEIGHT_OEM, warn );
#ifdef LIBIMP_DEBUG
                if ( warn ) {
                    ILog( EWarn, "There are weights smaller than %f!", MIN_WEIGHT_OEM );
                }
#endif
#else
                typedef TDataSet DatasetType;
                const TDataSet & dataset = inputDataset;            
#endif
                info.nAdjustedSamples = nDatasetSize;

                /// we need to know this because of MAP update 
                model.setDatasetSize( (unsigned int) nDatasetSize );

                ImStaticArray<TScalar,MAX_FITTED_LOBES> pdf;
//                int batchSize = clamp( cfg.batchSize, 1, nDatasetSize );
                        
                MixtureLikelihoodFunctional<TMixtureModel, DatasetType, TScalar, DIMENSION> loglkh( model );
                Float prevLkh, currLkh = -std::numeric_limits<Float>::max();
                size_t iterCounter = 0;
				const int lobesCount = model.nComponents() / DIMENSION;
				static const TScalar v1 = TScalar( 1.0f );

                do {
                    //std::random_shuffle( dataset.begin(), dataset.end() );
                    iterCounter ++;
                    prevLkh = currLkh;
                    for ( int i = 0; i < nDatasetSize; ++i ) {

						//////////////////////////////////////////////////////////////////////////
						/// E-step: update sufficient statistics
						//////////////////////////////////////////////////////////////////////////

						/// compute posterior probabilities                
						std::fill( pdf.begin(), pdf.begin() + lobesCount, TScalar(0.f) );
						TScalar denom				  = TScalar( 0.f );
						const TScalar & dataWeight    = dataset[ i ].getValue();            
						const TVector & x			  = dataset[ i ].getPoint();

						for ( int h = 0; h < lobesCount; ++h ) {
							pdf[ h ] = model.responsibility( h, x );
							denom += pdf[ h ];
						}
						denom = TScalar( TScalar::dot( denom, v1 ) );
						auto stats  = model.getStats();
						Float eta   = std::pow( stats.getCount() + 1.f, -cfg.alpha );
						/// update weight statistic (common for all components)
						stats.updateWeight( eta, dataWeight[ 0 ] );
						TScalar etaV( eta );

						for ( int h = 0; h < lobesCount; ++h ) {
							TScalar posterior = pdf[ h ] / denom;
							stats.updateMultiLobe( h, etaV, x, dataWeight, posterior );
						}            
						stats.increment();

						//////////////////////////////////////////////////////////////////////////
						/// M-step: update model
						//////////////////////////////////////////////////////////////////////////         
                
						//if ( model.getStats().getCount() + 1 > cfg.delayedUpdate ) {
						if ( model.getStats().getCount() % cfg.batchSize == 0 ) {
							for ( int h = 0; h < lobesCount; ++h ) {
								model.update( h, cfg.regularizationType );
							}
						}

						IMPORTANCE_ASSERT( model.check() );
					}

                    currLkh = loglkh.eval( dataset, nDatasetSize );

                } while ( std::abs( prevLkh / currLkh - 1.f ) >= TMixtureModel::getMinimalError() && iterCounter < cfg.nMaxEMIter );


                info.lkh = currLkh;
                info.nIter = (int) iterCounter;
#ifdef LIBIMP_STATS
                emStats.nEMIterations += iterCounter;
                emStats.nEMRuns++;
                model.getDistribution().info.nIteration = iterCounter;
#endif
        }
    }
}