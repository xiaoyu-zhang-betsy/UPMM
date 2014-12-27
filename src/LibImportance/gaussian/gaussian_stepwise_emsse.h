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
#pragma warning(disable:4512)

#include "../LibImportanceTypes.h"
#include "gaussiansse.h"

namespace Importance {
    namespace EM {

        template<typename TScalar, int DIMENSION, typename Distribution>
        class GaussianSufficientStats {
        public:
            GaussianSufficientStats( Distribution & dist ) : m_dist( dist ) {}

            // weight update does not depend on any component thus it is separate
            IMPORTANCE_INLINE void updateWeight( Float eta, Float weight ) {
                m_dist.statsPointWeight = ( 1.f - eta ) * m_dist.statsPointWeight  + eta * weight;
            }

            IMPORTANCE_INLINE void updateMultiLobe( int h, TScalar eta, const QuadVector2 & point, TScalar weight, TScalar posterior) {
                static const TScalar v1 = TScalar( 1.0f );
                GaussianMultiLobe< TScalar > & lobe = m_dist.lobes[ h ];
                
                TScalar & statsW					= lobe.statsW;
                TScalar & statsGamma_w        = lobe.statsGamma_w;
                
                TScalar a		= v1 - eta;
                TScalar b		= eta * weight * posterior;

                statsW                  = a * statsW                    + b;                
                lobe.x1                 = a * lobe.x1                   + b * point.x;
                lobe.x2                 = a * lobe.x2                   + b * point.y;
                lobe.x1x2               = a * lobe.x1x2                 + b * point.x * point.y;
                lobe.x1sqr              = a * lobe.x1sqr                + b * point.x * point.x;
                lobe.x2sqr              = a * lobe.x2sqr                + b * point.y * point.y;
                // we need to keep track of this statistic for covariance regularization (MAP solution)
                statsGamma_w              = (v1 - eta) * statsGamma_w  + eta * posterior;

                IMPORTANCE_ASSERT( statsW.isReal() );
                IMPORTANCE_ASSERT( lobe.x1.isReal() && lobe.x2.isReal() );
            }

            IMPORTANCE_INLINE unsigned int getCount() const { return m_dist.k; }

            IMPORTANCE_INLINE void increment() { m_dist.k ++; }
        private:
            Distribution & m_dist;            
        };

		template<class TScalar>
		inline void failCheck( const Bool4 & crit, const SPDMatrix2x2<TScalar> & trueVal, SPDMatrix2x2<TScalar> & out ) {
			out.m[ 0 ] = crit.blend( trueVal.m[ 0 ], out.m[ 0 ] );
			out.m[ 1 ] = crit.blend( trueVal.m[ 1 ], out.m[ 1 ] );
			out.m[ 2 ] = crit.blend( trueVal.m[ 2 ], out.m[ 2 ] );
		}
        
        template<class TScalar, int DIMENSION, class Distribution>
        class GaussianMixtureModel {
        public:
            GaussianMixtureModel( Distribution & dist ) : m_dist( dist ), m_datasetSize( std::numeric_limits<unsigned int>::max() ) {
				for ( int h = 0; h < m_dist.lobes.size(); ++h ) {
					weights[ h ] = m_dist.lobes[ h ].weights;
				}
			}

            IMPORTANCE_INLINE void update( int h, ERegularizationType regularizationType ) {
                GaussianMultiLobe< TScalar > & g = m_dist.lobes[ h ];                          
                static const TScalar v0 = TScalar( 0.0f );
				static const TScalar minFloat = TScalar( std::numeric_limits<Float>::min() );
                const TScalar dimV = TScalar( (float)m_dist.nComponents() );

                IMPORTANCE_ASSERT( (g.statsW > v0).allTrue() );            
				auto crit = g.statsW > minFloat;

                IMPORTANCE_ASSERT( m_dist.statsPointWeight > 0.f );
                g.mean.x			= crit.blend( g.x1 / g.statsW, g.mean.x );
				g.mean.y			= crit.blend( g.x2 / g.statsW, g.mean.y );
                
                TScalar m0 = g.x1sqr - 2*g.x1*g.mean.x                      + g.statsW * g.mean.x*g.mean.x;
                TScalar m1 = g.x1x2  - g.x1*g.mean.y - g.x2*g.mean.x        + g.statsW * g.mean.x*g.mean.y;
                TScalar m2 = g.x2sqr - 2*g.x2*g.mean.y                      + g.statsW * g.mean.y*g.mean.y;
                SPDMatrix2x2<TScalar> denormCov( m0, m1, m2 );                

                /// covariance                
                switch ( regularizationType ) {
                case EMAP: {                        
                    unsigned int n  = getStats().getCount();
                    n               = std::min( m_datasetSize, n );
                    TScalar nV		= TScalar((float)n);
                    const TScalar & a   = Distribution::PRIOR_A;                            
                    const TScalar & b	= Distribution::PRIOR_B;
                    const TScalar & v	= Distribution::PRIOR_V;                            

                    const SPDMatrix2x2<TScalar> covmat = (SPDMatrix2x2<TScalar>(b, v0, b) + nV * (denormCov / TScalar(m_dist.statsPointWeight))) 
                                                            / ( a + nV * (g.statsW / TScalar(m_dist.statsPointWeight)) );
                    failCheck( crit, covmat, g.cov ); 
                    g.weights	= crit.blend( ( ( g.statsW / TScalar( m_dist.statsPointWeight ) ) * nV + v ) / ( nV + dimV * v ), g.weights );
                }
                break;

                case ESatoIshy: {                        
                    IMPORTANCE_ASSERT( false );
                    g.weights			= crit.blend( g.statsW   /  TScalar( m_dist.statsPointWeight ), g.weights );
                    failCheck( crit, denormCov / g.statsW, g.cov );                    


                    const TScalar p = TScalar( 0.01f );
                    TScalar deltaSqr = g.cov.trace() * TScalar( 0.5f );
                    deltaSqr = TScalar( _mm_max_ps(_mm_set1_ps(1e-3f), deltaSqr.data) );
                    g.cov.m[ 0 ] += crit.blend( p * deltaSqr, v0 );
                    g.cov.m[ 2 ] += crit.blend( p * deltaSqr, v0 );
                }
                break;

                case ENone:
                    IMPORTANCE_ASSERT( false );
                    g.weights			= crit.blend( g.statsW   /  TScalar( m_dist.statsPointWeight ), g.weights );
                    failCheck( crit, denormCov / g.statsW, g.cov );
                break;

                default:
                    IMPORTANCE_ASSERT( false );
                }

                Bool4 degenerate = g.cov.determinant() < Distribution::MIN_DETERMINANT;

                // Restart
                if( !degenerate.allFalse() ) {
                    for ( int i = 0; i < DIMENSION; ++i )  {
                        if ( degenerate[ i ] ) {
                            g.cov.m[ 0 ][ i ] = 0.1f;
                            g.cov.m[ 1 ][ i ] = 0.0f;
                            g.cov.m[ 2 ][ i ] = 0.1f;
#ifdef LIBIMP_DEBUG
                            ILog( EWarn, "Gaussian online: Variance of component %d was reset due to overfitting!", h * DIMENSION + i );
#endif
                        }						
                    }
                }

                IMPORTANCE_ASSERT( g.mean.x.isReal() );

                g.init();
            }

            IMPORTANCE_INLINE int nComponents() const { return (int) m_dist.nComponents(); }

            IMPORTANCE_INLINE GaussianSufficientStats<TScalar, DIMENSION, Distribution> getStats() { return GaussianSufficientStats<TScalar, DIMENSION, Distribution>( m_dist ); }

            IMPORTANCE_INLINE void setDatasetSize( unsigned int val ) { m_datasetSize = val; }

            IMPORTANCE_INLINE TScalar responsibility( int h, const QuadVector2 & point ) const {
                IMPORTANCE_ASSERT( h < m_dist.lobes.size() );
                return m_dist.lobes[ h ].weights * m_dist.lobes[ h ].pdf( point );
            }

            static Float getMinimalError() {
                return 1e-3f;
            }

            /// debug 
            bool check() const { 
                Float error, weightsSum;
                return m_dist.check( error, weightsSum );
            }

            std::string toString() const { return m_dist.toString(); }
            Distribution & getDistribution() { return m_dist; }
        private:
            Distribution & m_dist;
			TScalar weights[ MAX_FITTED_LOBES / DIMENSION ];
            unsigned int m_datasetSize;
        };        

    }
}
#pragma warning(pop)