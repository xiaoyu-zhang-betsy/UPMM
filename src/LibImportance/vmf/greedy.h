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

#include "../shared/Config.h"
#include "../LibImportanceTypes.h"
#include "random.h"
#include "datasetrecord.h"
#include "Vmf.h"

#pragma warning(push)
#pragma warning(disable:4100)
#pragma warning(disable:4512)
#pragma warning(disable:4996)

namespace Importance {

    class HardAssignOrdering : public std::binary_function<size_t, size_t, bool> {
    public:
        inline HardAssignOrdering( const IStack<size_t> & hardAssignment ) : m_hs( hardAssignment ) {}
        inline bool operator()( size_t i1, size_t i2) const {
            return m_hs[ i1 ] < m_hs[ i2 ];
        }
    private:
        const IStack<size_t> & m_hs;
    };

    class DatasetVect {
    public:
        DatasetVect( size_t n ) : m_weight( 0.f ) { m_data.reserve( n ); }
        inline void add( const Particle & particle, const Vector3 & normal ) {
            m_weight += particle.weight;
            m_data.push( DatasetRecord( particle ) );
        }
        inline size_t getParticlesCount() const { return m_data.size(); }
        inline IStack<DatasetRecord> getDataset() { return m_data; }
        inline Float particlesTotalWeight() { return m_weight; }
        inline IStack<DatasetRecord>::iterator begin() {
            return m_data.begin();
        }
        inline IStack<DatasetRecord>::iterator end() {
            return m_data.end();
        }
    private:
        Float m_weight;
        IStack<DatasetRecord> m_data;
    };
    
    class GreedyToolSet {
    private:
        IMPORTANCE_INLINE void prepareCandidates( const VMF & motherComp, VMF * candidates, size_t m ) const {
            const Vector3 & dir = motherComp.m_mi;
            Float newKappa      = 1.25f * motherComp.m_kappa;
            
            Vector3 axis, t;
            coordinateSystem( dir, axis, t );
            Float clampedKappa = clamp( motherComp.m_kappa, 0.f, VMFM::MAX_KAPA );
            Float tmp = 1.f + std::log( 0.5f ) / clampedKappa;
            Float alpha = acos( tmp );
            Vector3 mi0 = Transform::rotation( axis, alpha ) * dir;
            
            candidates[ 0 ].m_mi    = mi0;
            candidates[ 0 ].m_kappa = newKappa;
            candidates[ 0 ].initPDFOnly();

            Float delta = ( 2.f * (Float) Importance::PI ) / m;
            Float angle = delta;
            for ( int i = 1; i < m; angle += delta, ++i ) {                
                candidates[ i ].m_mi    = Transform::rotation( dir, angle ) * mi0;
                candidates[ i ].m_kappa = newKappa;
                candidates[ i ].initPDFOnly();
            }
        }


        /* Implementation according to Verbeek, Vlasis: Efficient Greedy Learning of Gaussian Mixture Models; it can change passed number of candidates (can be only lowered) */
        IMPORTANCE_INLINE void prepareCandidates( const DatasetRecord * dataset, size_t n, VMF * candidates, size_t & m ) {
            IMPORTANCE_ASSERT( n > 1 );

            bool found = true;
            for ( size_t k = 0; found && k < m; ) {
                size_t index = std::min( (size_t) std::floor( n * nextFloat( m_rnd ) ), n - 1 );
                const Vector3 & xl = dataset[ index ].dir;
                index = std::min( (size_t) std::floor( n * nextFloat( m_rnd ) ), n - 1 );
                const Vector3 & xr = dataset[ index ].dir;

                Vector3 mil( 0.f ), mir( 0.f );
                Float nsl = 0.f, nsr = 0.f;
                for ( size_t i = 0; i < n; ++i ) {
                    const DatasetRecord & rec = dataset[ i ];
                    if ( dot( rec.dir, xl ) > dot( rec.dir, xr ) ) {
                        mil += rec.value * rec.dir;
                        nsl += rec.value;
                    } else {
                        mir += rec.value * rec.dir;
                        nsr += rec.value;
                    }
                }
                
                found = false;
                if ( k < m && nsl > 1e-7f ) {
                    candidates[ k++ ].computeParams( mil, nsl );
                    found = true;
                }
                if ( k < m && nsr > 1e-7f ) {
                    candidates[ k++ ].computeParams( mir, nsr );
                    found = true;
                }
            }
        }

        void prepareCandidatesKmeanspp( DatasetRecord * dataset, size_t nDatasetSize, Random & rnd, VMF * dists, Float * weights, size_t & m ) const {
            m = clamp( nDatasetSize / 10, (size_t) 1, (size_t) m );
            if ( nDatasetSize < 5 ) {
                m = 0;
                return;
            }            
            kmeanspp( dataset, nDatasetSize, rnd, dists, weights, m );
        }

        IMPORTANCE_INLINE void prepareCandidates( Float minimalConcentration, const DatasetRecord * dataset, size_t n, VMF * candidates, size_t m ) {
            for ( size_t i = 0; i < m; ++i ) {
                /// chose random subset of points
                size_t nSubs = std::min( (size_t) std::floor( n * nextFloat( m_rnd ) + 1 ), n );
                Vector3 mi( 0.f );
                Float ns = 0.f;

                // estimate parameters of a new candidate (component)
                for ( size_t j = 0; j < nSubs; ++j ) {
                    size_t index = std::min( (size_t) std::floor( n * nextFloat( m_rnd ) ), n - 1 );
                    const DatasetRecord & rec = dataset[ index ];
                    mi += rec.value * rec.dir;
                    ns += rec.value;
                }                
                candidates[ i ].computeParams( mi, ns );                                           
                
                if ( candidates[ i ].m_kappa < minimalConcentration ) {
                    candidates[ i ].m_kappa = minimalConcentration;
                    candidates[ i ].initPDFOnly();
                }
            }
        }

        IMPORTANCE_INLINE Float newCompLkhContrib( const VMF & newComp, const DatasetRecord & rec, Float a, Float ra ) const {
            if ( a < 1e-7f || ra < 1e-7f ) {
                return 0.f;
            }

            return rec.value * ra * ( std::log( a / ra ) + newComp.logAveraged( rec.dir ) );
        }
  
        IMPORTANCE_INLINE Float fixedMixtureLkhContrib( const DatasetRecord & rec, Float a, Float ra ) const {
            /// Note that ra depends on a. If a goes to 1 then ra goes to 1 as well and limit of the expression below goes to zero.                        
            if ( 1.f - a < 1e-7f || 1.f - ra < 1e-7f ) {
                return 0.f;
            }

            return ( 1.f - ra ) * ( rec.value * std::log(( 1.f - a )/(1.f - ra)) + rec.fa );
        }
    
        IMPORTANCE_INLINE bool partialEM( 
            const VMFM & fixedDist, VMF & newComp, Float & a, 
            const DatasetRecord * dataset, size_t n, 
            const DatasetRecord * zeroPosteriorDataset1, size_t zeroPosterSetSize1,
            const DatasetRecord * zeroPosteriorDataset2, size_t zeroPosterSetSize2,
            Float & lkh ) const {

            int dummy;
            return partialEM( fixedDist, newComp, a, dataset, n, 
                zeroPosteriorDataset1, zeroPosterSetSize1, 
                zeroPosteriorDataset2, zeroPosterSetSize2, 
                lkh, dummy );
        }
        /* 
         params:         
         fixedDist  - k component mixture which is not changed in the partial E-M
         newComp    - new k+1st component
         a          - weight of new compoennt in k+1 mixture
         dataset    - cellpoints (datapoints)
         n          - number of datapoints
         lkh        - log-likelikhood of the dataset under k-component mixture (fixedDist), output is log-likelihood under new k+1 component mixture
         iIter      - output param, information about number of iterations
        */
        IMPORTANCE_INLINE bool partialEM( 
            const VMFM & fixedDist, VMF & newComp, Float & a, 
            const DatasetRecord * dataset, size_t n, 
            const DatasetRecord * zeroPosteriorDataset1, size_t zeroPosterSetSize1,
            const DatasetRecord * zeroPosteriorDataset2, size_t zeroPosterSetSize2,
            Float & lkh, int & iIter ) const {

            IMPORTANCE_ASSERT( a >= 0.f && a <= 1.f );            
            lkh         = -std::numeric_limits<Float>::max();
            Float prevLkh;
            iIter = 0;            

            do {
                iIter++;
                prevLkh = lkh;
                lkh = 0.f;
                Float newA = 0.f;
                Float ns = 0.f;
                Vector3 newMi( 0.f );
                for ( size_t i = 0; i < n; ++i ) {
                    //E and M together in one loop
                    const DatasetRecord & rec = dataset[ i ];                    

                    // E
                    Float e     = Ff::importanceExp( rec.fa / rec.value );
                    Float tmp   = a * newComp.pdf( rec.dir );
                    Float ra    = ( tmp ) / ( ( 1.f - a ) * e + tmp );
                    IMPORTANCE_ASSERT( !isNaN( ra ) );

                    //M
                    ns      += rec.value;
                    newA    += rec.value * ra;
                    newMi   += rec.value * ra * rec.dir;

                    lkh += newCompLkhContrib( newComp, rec, a, ra );
                    lkh += fixedMixtureLkhContrib( rec, a, ra );     
                    IMPORTANCE_ASSERT( !isINF( lkh ) && !isNaN( lkh ) );               
                } 

                for ( size_t i = 0; i < zeroPosterSetSize1; ++i ) {
                    lkh += fixedMixtureLkhContrib( zeroPosteriorDataset1[ i ], a, 0.f );
                    IMPORTANCE_ASSERT( !isINF( lkh ) && !isNaN( lkh ) );
                }
                for ( size_t i = 0; i < zeroPosterSetSize2; ++i ) {
                    lkh += fixedMixtureLkhContrib( zeroPosteriorDataset2[ i ], a, 0.f );
                    IMPORTANCE_ASSERT( !isINF( lkh ) && !isNaN( lkh ) );
                }                
                
                newComp.computeParams( newMi, newA );
                a = newA / ns;

            } while ( std::abs( prevLkh / lkh - 1.f ) >= VMFM::MINIMAL_ERROR && lkh > prevLkh && iIter < 40 );
            //ILog( EInfo, "a = %f", a );
            return a > 1e-7f;
        }

        struct Stats {
            Stats() : mi( 0.f ), ns( 0.f ), n( 0 ) {}
            Vector3 mi;
            Float   ns;
            size_t  n;
        };

    public:
        static IMPORTANCE_INLINE void kmeanspp( const DatasetRecord * dataset, size_t nDatasetSize, Random & rnd, VMF * dists, Float * weights, size_t m ) {
            const int MAX = 300;

            if ( m == 0 ) {
                return;
            }

            IMPORTANCE_ASSERT( m < MAX );
            IMPORTANCE_ASSERT( nDatasetSize > 0 );
            IMPORTANCE_ASSERT( nDatasetSize > m );


            ImStaticArray<Stats,MAX> stats;
            ImStaticArray<Vector3,MAX> centers;           
            Float * pdf  = (Float*) alloca( ( nDatasetSize + 1 ) * sizeof( Float ) ); 

            

            /// chose first centroid randomly
            size_t index = std::min( (size_t) std::floor( nDatasetSize * nextFloat( rnd ) ), nDatasetSize - 1 );
            centers[ 0 ] = dataset[ index ].dir;
            
            /// add new centroids up to number m
            for ( int i = 1; i < m; ++i ) {
                Float sum = 0.f;

                /// search for candidates over all datapoints (we set probability to each of them)
                for ( int j = 0; j < (int) nDatasetSize; ++j ) {                    
                    pdf[ j ] = -1.f;
                    // check distance to all selected centroids and pick the closest one
                    for ( int k = 0; k < i; ++k ) {
                        Float dotProductDist = dot( centers[ k ], dataset[ j ].dir );
                        if ( dotProductDist > pdf[ j ] ) {
                            pdf[ j ] = dotProductDist;
                        }
                    }
                    // transform cosine distance so that probability of selecting a point with closest neighbour is the least
                    pdf[ j ] = std::max( pdf[ j ] * -1.f + 1.f, 0.f );
                    sum += pdf[ j ];
                    IMPORTANCE_ASSERT( pdf[ j ] >= 0.f && pdf[ j ] < 2.0f );
                }

                /// pick new centroid according to pdf constructed just above
                Float u             = nextFloat(rnd) * sum;            
                Float cdfVal        = pdf[ 0 ];            
                pdf[ nDatasetSize ] = 0.f;
                size_t l;
                for ( l = 0; l < nDatasetSize && cdfVal <= u; ++l ) {
                    cdfVal += pdf[ l + 1 ];
                }
                IMPORTANCE_ASSERT( l < nDatasetSize && cdfVal > u );
                centers[ i ] = dataset[ l ].dir;
            }

            // compute parameters of distributions according to factorization determined by centroids
            Float totalWeight = 0.f;
            for ( int j = 0; j < (int) nDatasetSize; ++j ) {
                const DatasetRecord & rec   = dataset[ j ];
                Float minDotDist            = -1.f;
                int factorIndex             = 0;
                /// find out the nearest centroid to current data point with index j
                for ( int i = 0; i < m; ++i ) {                    
                    Float dotDist = dot( rec.dir, centers[ i ] );
                    if ( minDotDist < dotDist ) {
                        minDotDist     = dotDist;
                        factorIndex    = i;
                    }
                }
                stats[ factorIndex ].n ++;
                stats[ factorIndex ].mi += rec.value * rec.dir;
                stats[ factorIndex ].ns += rec.value;                
                totalWeight             += rec.value;
            }

            for ( int i = 0; i < m; ++i ) {
                const Stats & s = stats[ i ];
                if ( weights != NULL ) {
                    weights[ i ] = s.ns / totalWeight;
                }
                IMPORTANCE_ASSERT( s.ns > 0.f ); // this should never happen as every centroid is closest to itself
                dists[ i ].computeParams( s.mi, s.ns );
                
                /// be somehow defensive to over-fitting
                if ( s.n < 4 && dists[ i ].m_kappa >= VMFM::MAX_KAPA ) {
                    dists[ i ].m_kappa = 10.f;
                    dists[ i ].initPDFOnly();
                }
            }
        }


        IMPORTANCE_INLINE static Float AIC( unsigned int n, unsigned int k, Float logLkh ) {
            return 2.f * ( k - logLkh );
        }

        IMPORTANCE_INLINE static Float AICc( unsigned int n, unsigned int k, Float logLkh ) {
            Float tmp = (Float) n-k-1;
            if ( tmp == 0.f ) {
                tmp = 1e-6f;
            }
            return AIC( n, k, logLkh ) + 2.f*k*(k+1)/tmp;
        }

        bool insertFirstComponent( VMFM & dist,
            IStack<DatasetRecord> & dataset,
            Float & lkh,
            size_t m = VMFM::M_CANDIDATES_GREEDY ) {
                IMPORTANCE_ASSERT( dist.nComponents() == 0 );
                IMPORTANCE_ASSERT( dataset.size() > 0 );

                if ( dist.c.size() == 0 ) {
                    dist.c.push( VMF() );
                    dist.w.push( 0.f );
                }

                VMF * candidates = static_cast<VMF*>( 
                    _alloca( sizeof( VMF ) * m ) );

                prepareCandidates( 0.f, &dataset[ 0 ], dataset.size(), candidates, m ); 

                Float maxLkh    = -std::numeric_limits<Float>::max();
                bool found      = false;                
                for ( size_t iCand = 0; iCand < m; ++iCand ) {
                    Float tmpLkh = -std::numeric_limits<Float>::max();
                    Float tmpA   = 1.f;
                    bool res = partialEM( dist, candidates[ iCand ], tmpA, 
                        &dataset[ 0 ], dataset.size(), NULL, 0, NULL, 0, tmpLkh );
                    // pick the best candidate according to log-likelihood as a new initial component
                    if ( res && tmpLkh > maxLkh ) {
                        dist.c[ 0 ] = candidates[ iCand ];
                        dist.w[ 0 ] = 1.f;
                        maxLkh      = tmpLkh;
                        found       = true;
                    }
                }

                if ( found ) {
                    dist.m_nComponents++;
                    lkh = maxLkh;
                }

                return found;
        }



        bool insertNewComponent( VMFM & dist, 
            IStack<DatasetRecord> & dataset, 
            const IStack<size_t> & hardAssignVect,
            Float & lkh, 
            size_t m = VMFM::M_CANDIDATES_GREEDY ) {
                        
                IMPORTANCE_ASSERT( dataset.size() == hardAssignVect.size() );                
                
                size_t k = dist.nComponents();
                IMPORTANCE_ASSERT( k  > 0 );
                if ( k >= dist.c.size() ) {
                    dist.c.push( VMF() );
                    dist.w.push( 0.f );
                }

                // factorization of points according to their nearest components
                IStack<size_t> factor( dataset.size() );
                /// work with indices so that sorting does not need to move around all data 
                for ( size_t i = 0; i < factor.size(); ++i ) { factor[ i ] = i; }
                /// factorize according to hard assignment
                std::sort( factor.begin(), factor.end(), HardAssignOrdering( hardAssignVect ) );

                
                
                /// allocate memory
                DatasetRecord * factorDataset = static_cast<DatasetRecord*>( 
                    _alloca( sizeof( DatasetRecord ) * dataset.size() ) );
                VMF * candidates = static_cast<VMF*>( 
                    _alloca( sizeof( VMF ) * m ) );
                


                /// prepare data points so that they are organized in one array in coherent bunches - every bunch corresponds to factor according to component h
                /// size of every bunch is stored in factorSize array
                Importance::ImStaticArray<int,Importance::MAX_FITTED_LOBES> factorSize;                
                int h        = -1;
                int n        = -1;
                int nFactors = 0;
                for ( size_t i = 0; i < factor.size(); ++i ) {                    
                    if ( (int) hardAssignVect[ factor[ i ] ] > h ) {
                        h = (int) hardAssignVect[ factor[ i ] ];                        
                        if ( n > 0 ) {                            
                            factorSize[ nFactors++ ] = n;
                        }
                        n = 0;
                    }
                    factorDataset[ i ] = dataset[ factor[ i ] ];
                    n++;
                }
                if ( n > 0 ) {                    
                    factorSize[ nFactors++ ] = n;
                }
                //IMPORTANCE_ASSERT( k == nFactors );


                Float aicc                  = AICc( (unsigned int) dataset.size(), dist.nParameters(), lkh );
                Float maxLkh                = lkh;
                bool found                  = false;   
                DatasetRecord * setBegin    = factorDataset;
                size_t zeroPosterSetSize1   = 0;
                for ( int q = 0; q < nFactors; q++ ) {

                    size_t setSize              = (size_t) factorSize[ q ];
                    size_t tmpM                 = m;
                    /// initialize set of m new candidates based on component h
                    //prepareCandidates( /*dist.c[ q ].m_kappa * 1.25f*/0.f, setBegin, setSize, candidates, tmpM );   
                    //prepareCandidates( dist.c[ q ], candidates, tmpM );
                    //prepareCandidates( setBegin, setSize, candidates, tmpM );
                    prepareCandidatesKmeanspp( setBegin, setSize, m_rnd, candidates, NULL, tmpM );
                    
                    /// weight of new candidate component in new k+1 mixture
                    Float a = dist.w[ q ] * 0.5f;    
                    /// apply partial E-M to all m candidates and pick the best one with respect to maximal log-likelihood of two component mixture
                    for ( size_t iCand = 0; iCand < tmpM; ++iCand ) {
                        Float tmpLkh    = lkh;
                        Float tmpA      = a;                      
                        int nIterations;
                        bool res = partialEM( dist, candidates[ iCand ], tmpA, 
                                                setBegin, setSize, 
                                                factorDataset, zeroPosterSetSize1,
                                                setBegin + setSize, factor.size() - zeroPosterSetSize1 - setSize,
                                                tmpLkh, nIterations );                       
                        //ILog( EInfo, "partial E-M iterations: %d", nIterations );
                        // pick the best candidate according to log-likelihood as a new k+1st component
                        if ( res && tmpLkh > maxLkh && 1.f - tmpA > 1e-6f ) {
                            dist.c[ k ] = candidates[ iCand ];
                            dist.w[ k ] = tmpA;
                            maxLkh      = tmpLkh;
                            found       = true;
                            //ILog( EInfo, "Accepted" );
                        }
                    }

                    zeroPosterSetSize1  += setSize;
                    setBegin            += setSize;
                }
                
                /// if there was a gain in inserting a new component and its weight does not eliminates the rest of a mixture
                found = found && ( 1.f - dist.w[ k ] ) > 1e-7f;
                if ( found ) {
                    dist.m_nComponents++;
                    lkh = maxLkh;
                    /// check AICc
                    Float tmp = AICc( (unsigned int) dataset.size(), dist.nParameters(), lkh );
                    if ( tmp < aicc && std::abs( tmp/aicc - 1.f ) > 1e-3f ) {                        
                        for ( size_t j = 0; j < k; ++j ) { 
                            dist.w[ j ] *= ( 1.f - dist.w[ k ] ); 
                        }
                    } else {
                        dist.m_nComponents--;
                        found = false;
                    }
                }

                return found;
        }

    private:
        Random m_rnd;
    };

}

#pragma warning(pop)