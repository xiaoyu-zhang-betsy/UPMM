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
#define ADJUST_DATASET 0

namespace Importance {  

    namespace OnlineEM {
        const int MAX_DATASETSIZE = 4 * MAX_KNN_PARTICLES;

        template<typename TVector>
        struct Sample {
            IMPORTANCE_INLINE const TVector & getPoint() const { return m_point; }
            IMPORTANCE_INLINE Float getValue() const { return m_weight; }
            IMPORTANCE_INLINE void setValue( const Float & w ) { m_weight = w; }
            //IMPORTANCE_INLINE Float getInverseDistance() const { return m_inverseDistance; }
            IMPORTANCE_INLINE Float getDistance() const { return m_inverseDistance; }

            TVector m_point;
            Float m_weight, m_inverseDistance;
        };

        template<typename TDataSet, typename TVector>
        IMPORTANCE_INLINE int adjustDataset( const TDataSet & dataset, int n, ImStaticArray<OnlineEM::Sample<TVector>,MAX_DATASETSIZE> & samples ) {
            IMPORTANCE_ASSERT( n > 0 );

            const Float minw = 0.f;

            // compute average value
            Float maxw       = 0.f;
            for ( int i = 0; i < n; ++ i ) {
                maxw += dataset[ i ].getValue();
            }
            maxw /= n;

            int j = 0;
            for ( int i = 0; i < n; ++ i ) {
                Float q = dataset[ i ].getValue();
                if ( q >= minw ) {                                                                
                    while ( q > maxw ) {
                        IMPORTANCE_ASSERT( j < MAX_DATASETSIZE );
                        q -= maxw;
                        samples[ j ].m_point           = dataset[ i ].getPoint();
                        samples[ j ].m_inverseDistance = dataset[ i ].getDistance();//getInverseDistance();
                        samples[ j ].setValue( maxw );
                        j ++;                        
                    }
                    if ( q >= minw ) {
                        IMPORTANCE_ASSERT( j < MAX_DATASETSIZE );
                        samples[ j ].m_point           = dataset[ i ].getPoint();
                        samples[ j ].m_inverseDistance = dataset[ i ].getDistance();//getInverseDistance();
                        samples[ j ].setValue( q );
                        j ++;
                    }
                }
            }

            std::random_shuffle( samples.begin(), samples.begin() + j );

            return j;
        }

        template<class TMixtureModel, class TDataSet, class TVector>
        void stepwiseEM( 
            const TDataSet & inputDataset,
            int nDatasetSize,                        
            TMixtureModel & model, 
            const TVector &,
            const StepwiseEMConfig & cfg = StepwiseEMConfig() ) { 
                IMPORTANCE_ASSERT( nDatasetSize > 0 );
                IMPORTANCE_ASSERT( cfg.alpha > 0.5f && cfg.alpha <= 1.f );            
           
#if ADJUST_DATASET
                ImStaticArray<Sample<TVector>,MAX_DATASETSIZE> dataset;
                nDatasetSize = adjustDataset<TDataSet,TVector>( inputDataset, nDatasetSize, dataset );
#else
                const TDataSet & dataset = inputDataset;
#endif

                ImStaticArray<Float,MAX_FITTED_LOBES> pdf;
                int batchSize = clamp( cfg.batchSize, 1, nDatasetSize );

                int i = 0;
                while ( i < nDatasetSize ) {

                    //////////////////////////////////////////////////////////////////////////
                    /// E-step: update sufficient statistics
                    //////////////////////////////////////////////////////////////////////////

                    /// Process the whole batch at once
                    for ( int j = 0; j < batchSize && i < nDatasetSize; ++j, ++i ) {
                        /// compute posterior probabilities                
                        std::fill( pdf.begin(), pdf.begin() + model.nComponents(), 0.f );
                        Float denom                 = 0.f;
                        const Float & dataWeight    = dataset[ i ].getValue();            
                        const TVector & x           = dataset[ i ].getPoint();
                        Float invDist               = dataset[ i ].getDistance();//getInverseDistance();

                        for ( int h = 0; h < model.nComponents(); ++h ) {
                            pdf[ h ] = model.responsibility( h, x );
                            denom += pdf[ h ];
                        }

                        auto stats  = model.getStats();
                        Float eta   = std::pow( stats.getCount() + 1.f, -cfg.alpha );
                        /// update weight statistic (common for all components)
                        stats.updateWeight( eta, dataWeight );

                        for ( int h = 0; h < model.nComponents(); ++h ) {
                            Float posterior = pdf[ h ] / denom;
                            stats.updateComponent( h, eta, x, dataWeight, posterior, invDist );
                        }            
                        stats.increment();
                    }

                    //////////////////////////////////////////////////////////////////////////
                    /// M-step: update model
                    //////////////////////////////////////////////////////////////////////////         
                
                    if ( model.getStats().getCount() + 1 > cfg.delayedUpdate ) {
                        for ( int h = 0; h < model.nComponents(); ++h ) {
                            model.update( h, cfg.regularizationType );
                        }
                    }

                    IMPORTANCE_ASSERT( model.check() );
                }
        }



        /// General way to compute log-likelihood of a given mixture (Gaussian, vMFm, ...)
        template<class TMixtureModel, class TDataSet>
        class MixtureLikelihoodFunctional {
            void operator=(const MixtureLikelihoodFunctional&);
        public:
            MixtureLikelihoodFunctional( const TMixtureModel & model ) 
                : m_model( model ) {}

            /// Compute log-likelihood for given dataset
            Float eval( const TDataSet & dataset, int nDatasetSize ) const {
                Float res = 0;
                for ( int i = 0; i < nDatasetSize; ++i ) {

                    Float val = 0.f;
                    for ( int h = 0; h < m_model.nComponents(); ++h ) {                 
                        val += m_model.responsibility( h, dataset[ i ].getPoint() );
                    }

                    res += dataset[ i ].getValue() * Ff::log( val );
                }

                return res;
            }

        private:
            const TMixtureModel & m_model;        
        };

        template<class TMixtureModel, class TDataSet, class TVector>
        void stepwiseStaticEM( 
            const TDataSet & inputDataset,
            int nDatasetSize,                        
            TMixtureModel & model, 
            const TVector &,
            EMStats & stats,
            Float & likelihood,
            const StepwiseEMConfig & cfg = StepwiseEMConfig() ) { 
                IMPORTANCE_ASSERT( nDatasetSize > 0 );
                IMPORTANCE_ASSERT( cfg.alpha > 0.5f && cfg.alpha <= 1.f );         

#if ADJUST_DATASET
                typedef ImStaticArray<Sample<TVector>,MAX_DATASETSIZE> DatasetType;
                DatasetType dataset;
                nDatasetSize = adjustDataset<TDataSet,TVector>( inputDataset, nDatasetSize, dataset );
#else
                typedef TDataSet DatasetType;
                const TDataSet & dataset = inputDataset;            
#endif

                /// we need to know this because of MAP update 
                model.setDatasetSize( (unsigned int) nDatasetSize );

                ImStaticArray<Float,MAX_FITTED_LOBES> pdf;
                int batchSize = clamp( cfg.batchSize, 1, nDatasetSize );
                        
                MixtureLikelihoodFunctional<TMixtureModel, DatasetType> loglkh( model );
                Float prevLkh, currLkh = -std::numeric_limits<Float>::max();
                size_t iterCounter = 0;

                do {
                    //std::random_shuffle( dataset.begin(), dataset.end() );
                    iterCounter ++;
                    prevLkh = currLkh;
                    int i   = 0;
                    while ( i < nDatasetSize ) {

                        //////////////////////////////////////////////////////////////////////////
                        /// E-step: update sufficient statistics
                        //////////////////////////////////////////////////////////////////////////

                        /// Process the whole batch at once
                        for ( int j = 0; j < batchSize && i < nDatasetSize; ++j, ++i ) {
                            /// compute posterior probabilities                
                            std::fill( pdf.begin(), pdf.begin() + model.nComponents(), 0.f );
                            Float denom                 = 0.f;
                            const Float & dataWeight    = dataset[ i ].getValue();            
                            const TVector & x           = dataset[ i ].getPoint();
                            Float invDist               = dataset[ i ].getDistance();//getInverseDistance();

                            for ( int h = 0; h < model.nComponents(); ++h ) {
                                pdf[ h ] = model.responsibility( h, x );
                                denom += pdf[ h ];
                            }

                            auto stats  = model.getStats();
                            Float eta   = std::pow( stats.getCount() + 1.f, -cfg.alpha );
                            /// update weight statistic (common for all components)
                            stats.updateWeight( eta, dataWeight );

                            for ( int h = 0; h < model.nComponents(); ++h ) {
                                Float posterior = pdf[ h ] / denom;
                                stats.updateComponent( h, eta, x, dataWeight, posterior, invDist );
                            }            
                            stats.increment();
                        }

                        //////////////////////////////////////////////////////////////////////////
                        /// M-step: update model
                        //////////////////////////////////////////////////////////////////////////         

                        if ( model.getStats().getCount() + 1 > cfg.delayedUpdate ) {
                            for ( int h = 0; h < model.nComponents(); ++h ) {
                                model.update( h, cfg.regularizationType );
                            }
                        }

                        IMPORTANCE_ASSERT( model.check() );
                    }

                    currLkh = loglkh.eval( dataset, nDatasetSize );

                } while ( std::abs( prevLkh / currLkh - 1.f ) >= GaussianMixture::MINIMAL_ERROR && iterCounter < cfg.nMaxEMIter );

                likelihood = currLkh;
                stats.nEMIterations += iterCounter;
                stats.nEMRuns++;
        }
    }
}