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


#ifndef ___VMF___
#define ___VMF___

#include "random.h" 
#include "../LibImportanceTypes.h"
#include "../distribution/DefaultDistributionModel.h"
#include "../shared/Utils.h"

#ifdef TWOPI
#undef TWOPI
#endif

#pragma warning(push)
#pragma warning(disable:4265)
#pragma warning(disable:4127)
#pragma warning(disable:4100)

#include "../shared/Vector3.h"
#include "../shared/frame.h"
#include "../shared/StaticArray.h"
#include "../shared/PointKdTreeShared.h"
#include "../shared/Bitmap.h"


#include <vector>
#include <sstream>

namespace Importance {

    class VMFM;

    class VMF {
    public:
        explicit VMF(const Vector3 & mi = Vector3(0.0f), Float kapa = 0.0f)
            : m_initialized( false ), m_pdfInitialized( false ) {
            m_mi = mi;
            m_kappa = kapa;

            m_harmonicDistance = 0.f;            
            statsW             = 0.f;
            statsMi            = Vector3( 0.f );            
        }        

        /** Returns probability of sampling given direction. */
        IMPORTANCE_INLINE Float pdf( const Vector3 & dir ) const {
            IMPORTANCE_ASSERT( m_pdfInitialized );
            if (m_kappa == 0.f) {
                return UNIFORM_PDF;
            }

            const Float cosAngle = dot(dir, m_mi);        
            const Float exponent = m_kappa * ( cosAngle - 1.0f );
            const Float eTerm = Ff::importanceExp(exponent);
            const Float res = m_invDenom * m_kappa * eTerm;
            return res;
        }

        /** Sample a direction given random sample in unit square. */
        IMPORTANCE_INLINE Vector3 sampleDirection( Float u1, Float u2 ) const {
            IMPORTANCE_ASSERT( m_initialized );    
            if (m_kappa == 0.0) {
                return m_coordSys.toWorld( squareToSphere( u1, u2 ) );
            }

            Float b = m_invKapa * std::log( m_expMinusTwoKapa * (1.0f - u1) + u1 );
            IMPORTANCE_ASSERT( b < 1e-6 );
            Float q = 1.0f + std::max( b, Float(-2.0f));

            Float c = std::sqrt(1.0f - q * q);
            Float angle = TWOPI * u2;      
    
            return m_coordSys.toWorld(
                Vector3( c * cos( angle ), 
                        c * sin( angle ),
                        q
                ));
        }


        /** Compute averaged log function over all photons in the cell. */
        IMPORTANCE_INLINE Float logAveraged( const Vector3 dir ) const {
            IMPORTANCE_ASSERT( m_pdfInitialized );
            if (m_kappa == 0.0) {
                /* 1/(4*PI) is a limit at kappa = 0 of the computed expression. */
                return UNIFORM_PDF;
            }

            /* Beware that the vector dir is expected to be the average over the cell of directions and does not 
               need to be normalized. */
            Float cosAngle = dot(dir, m_mi);            
            const Float exponent = m_kappa * ( cosAngle - 1.0f );
            return std::log( m_invDenom * m_kappa ) + exponent;
        }


        /** Prepares the distribution for both sampling and computing PDF. */
        IMPORTANCE_INLINE void init() {                
            //you need this for computing the pdf (and part of it for sampling)
            initPDFOnly();
            //the rest of the initialization for sampling
            initForSampling();

            m_initialized = true;
        }	

    
        /** Prepares function only for computing PDF. It saves unnecessary constant initializations 
         *  which are needed only when we need to sample from the distribution. During EM the PDF
         *  is computed extensively but no sampling is needed.
         */
        IMPORTANCE_INLINE void initPDFOnly() {
            //you need this for both sampling and computing pdf
            m_expMinusTwoKapa = std::pow( EXP_MINUSTWO, m_kappa );

            //you need this for computing pdf only
            m_invDenom = 1.0f / (TWOPI * ( 1.0f - m_expMinusTwoKapa ));

            m_pdfInitialized = true;
        }


        /**
         * Initialize the rest of constants which are needed for sampling from the distribution. 
         * Prerequisite: initPDFOnly() must have been already called.
         */
        IMPORTANCE_INLINE void initForSampling() {
            IMPORTANCE_ASSERT( m_pdfInitialized );
            //you need this for sampling from this distribution only
            m_invKapa = 1.0f / (m_kappa > 0.0f ? m_kappa : 1.0f);
            m_coordSys = Frame(m_mi);

            m_initialized = true;
        }

        IMPORTANCE_INLINE Float normalizationConstant() const {
            if(m_kappa == 0) {
                return 1.f/(4*Importance::PI);
            } else {
                return m_invDenom * m_kappa;
            }
        }

        /// TODO: use this everywhere
        IMPORTANCE_INLINE void computeParams( const Vector3 & unnormalizedMi, Float weight, Float maxKappa = std::numeric_limits<Float>::max() ) {
            IMPORTANCE_ASSERT( weight > std::numeric_limits<Float>::min() );

            Float length    = unnormalizedMi.length();
            length          = std::max( length, std::numeric_limits<Float>::min() );
            m_mi            = unnormalizedMi / length;
            Float r         = length / weight;
#ifdef IMPORTANCE_SINGLE_PRECISION
            r                       = std::min( r, Float(0.99999997f) );
#else
            /// TODO: check this for double - hasn't been used yet
            IMPORTANCE_ASSERT( false );
            r                       = std::min( r, Float(0.999999999997f) );
#endif
            Float rSqr    = r * r;
            m_kappa = r * ( 3.0f - rSqr ) / ( 1.0f - rSqr );


            /// Checks 

            if ( !(m_kappa >= 0.f) || Importance::isNaN( m_kappa ) || Importance::isINF( m_kappa ) ) {
                ILog( EInfo, "length: %f; mi: %s; r: %f, weight: %f", length, m_mi.toString().c_str(), r, weight );
            }

            /// clamping
            if ( m_kappa > maxKappa ) {
                m_kappa = maxKappa;
            }
    
            IMPORTANCE_ASSERT( m_kappa >= 0.f );
            IMPORTANCE_ASSERT( !Importance::isNaN( m_kappa ) 
                            && !Importance::isINF( m_kappa ) ); 

            IMPORTANCE_ASSERT( m_mi.isNormalized() );

            initPDFOnly();
        }

        std::string toString() const {
            std::stringstream str;
            str << "Kapa: " << m_kappa << " Mi: " << m_mi.toString() << "HarmDist: " << m_harmonicDistance;
            return str.str();
        }

    public:
        Vector3 m_mi;
        Float m_kappa;
        Frame m_coordSys;        
        static const Float UNIFORM_PDF;
        static const Float TWOPI;
        static const Float EXP_MINUSTWO;
        //precomputed values dependent on kapa
        Float m_expMinusTwoKapa;
        Float m_invDenom;
        Float m_invKapa;
        bool m_initialized;
        bool m_pdfInitialized;
        Float m_harmonicDistance;
        /// Statistics for stepwise EM
        Float statsW;        
        Vector3 statsMi;
    };

    class VMFM : public DefaultDistributionModel {
    private:
        static const Float EPSILON;        
    public:    
        static const Float MINIMAL_STEP;
        static const Float MINIMAL_ERROR;
        static const Float MAX_KAPA;
        static const size_t MINIMAL_DATASIZE = 5;
        static const size_t M_CANDIDATES_GREEDY = 20;

    public:

        typedef VMF LobeType;


        VMFM() : m_initialized( false ) {
            m_nComponents    = 0;
            k                = 0;
        }


        explicit VMFM( size_t nComponents ) : m_initialized( false ) {
            m_initialized               = false;
            m_nComponents               = nComponents;
            k                           = 0;

            c.resize( m_nComponents, VMF() );
            w.resize( m_nComponents, 0.0f );
            cdf.resize( m_nComponents + 1, 0.0f );
        }

        
        IMPORTANCE_INLINE static float klDistance(const VMFM& a, const VMFM& b) {
            IMPORTANCE_ASSERT(false);
            return NAN;
        }
        IMPORTANCE_INLINE static float lobeDistance(const VMFM& a, const int indexA, const VMFM& b, const int indexB) {
            const float meanSqr = dot(a.c[indexA].m_mi, b.c[indexB].m_mi);
            const float wDiff = a.w[indexA] - b.w[indexB];
            const float matrixDiff = a.c[indexA].m_kappa - b.c[indexB].m_kappa;
            return Ff::sqrtFast(meanSqr + wDiff*wDiff + matrixDiff*matrixDiff);
        }
        IMPORTANCE_INLINE static LobeType lobeLerp(const VMFM& a, const int indexA, const VMFM& b, const int indexB, const float amountB, float& outW) {
            VMF result(lerp(a.c[indexA].m_mi, b.c[indexB].m_mi, amountB).getNormalized(),
                       lerp(a.c[indexA].m_kappa, b.c[indexB].m_kappa, amountB));
            outW = lerp(a.w[indexA], b.w[indexB], amountB);
            return result;
        }


        IMPORTANCE_INLINE Float pdf( const Vector3 & dir ) const {
            IMPORTANCE_ASSERT( m_initialized );
            Float res = 0.0;
            for (size_t i = 0; i < nComponents(); ++i ) {
                res += w[ i ] * c[ i ].pdf(dir);
            }

            return res;
        }
        IMPORTANCE_INLINE Importance::Vector3 sampleDirection( const Importance::Vector2 & sample ) const {
            IMPORTANCE_ASSERT( m_initialized );
            size_t index = getSampledComponentIndex( sample.x );
            IMPORTANCE_ASSERT( index < nComponents() );
            /// reuse sample here
            Float reusedSample = ( sample.x - cdf[ index ] ) 
                / ( cdf[ index + 1 ] - cdf[ index ] );
            IMPORTANCE_ASSERT( reusedSample >= 0.f && reusedSample <= 1.f );
            Vector3 res = c[ index ].sampleDirection( reusedSample, sample.y );
            return Importance::Vector3( res.x, res.y, res.z );
        }

        virtual const EmInfo * getInfo() const { return &m_info; }

        void fillBitmap(ImBitmap<float>& output) const {
            for(int x = 0; x < output.getWidth(); x++) {   
                for(int y = 0; y < output.getHeight(); y++) {
                    output(x, y) = this->pdf(mapDirection(x, y, output.getWidth(), output.getHeight()));
                }
            }
        }

        IMPORTANCE_INLINE void init() {
            cdf.resize( nComponents() + 1, 0.f );
            //compute cdf        
            Float sum = 0.0f;
            /// cdf is postponed by one index for comfy samples reusing
            cdf[ 0 ] = 0.f;
            for( size_t h = 1; h < nComponents() + 1; ++h ) {            
                sum += w[ h - 1 ];
                cdf[ h ] = sum;
                //finish initialization of each component
                c[ h - 1 ].initForSampling();
            }
            IMPORTANCE_ASSERT( std::abs(1.0f - cdf[ nComponents() ]) < EPSILON );
            cdf[ nComponents() ] = 1.0f;        
            m_initialized = true;
        }   

        void randomInitHemi( Importance::Random & rnd, const Vector3 & nrm );


        /** For debugging purposes - you can omit invocation of this method later for optimalization */
        IMPORTANCE_INLINE void invalidate() {
            m_initialized = false;
            for ( size_t h = 0; h < nComponents(); ++h ) {
                c[ h ].m_initialized = false;
                c[ h ].m_pdfInitialized = false;
            }
        }

        /* Check weights - have to sum to one; Check kappas - must be in specified interval */
        IMPORTANCE_INLINE bool check( Float & error, bool & isKappasValid, Float & weightsSum ) {        
            weightsSum      = 0.0f;   
            isKappasValid   = true;
            for ( size_t h = 0; h < nComponents(); ++h ) {
                IMPORTANCE_ASSERT(c[h].m_mi.isReal());
                weightsSum += w[ h ];
                if ( c[ h ].m_kappa < 0.0f /*|| c[ h ].m_kappa > MAX_KAPA*/ ) {
                    isKappasValid = false;
                }
            }

            error = std::abs( 1.0f - weightsSum );
            return isKappasValid && error < EPSILON; ;
        }

        /** Returns a number of components in the mixture. */
        IMPORTANCE_INLINE size_t nComponents() const {
            return m_nComponents;
        }

        IMPORTANCE_INLINE void initAsUninformed() {
            initAsUninformed(Vector3(0, 0, 1), 1.f);
        }

        IMPORTANCE_INLINE void initAsUninformed( const Vector3 & normal, Float harmonicDistance ) {
            IMPORTANCE_ASSERT( nComponents() == 0 );
            w.push( 1.f );
            VMF vmf;
            vmf.m_kappa             = 10.f;
            vmf.m_mi                = normal;
            vmf.m_harmonicDistance  = harmonicDistance;
            vmf.initPDFOnly();
            c.push( vmf );
            m_nComponents           = 1;
        }

        IMPORTANCE_INLINE void reserve( size_t nComponents ) {
            c.reserve( nComponents );
            w.reserve( nComponents );
            cdf.reserve( nComponents + 1 );
        }

        /* For use of AIC (Akaike Information Criterion) */
        IMPORTANCE_INLINE unsigned int nParameters() const {
            return (unsigned int) nComponents() * 10;
        }

        //void testCDF();

        virtual std::string toString() const {
            std::stringstream str;
            for ( size_t h = 0; h < nComponents(); ++h ) {
                str << "\tC[ " << h << " ]: " << "w: "<< w[ h ] << " " << c[ h ].toString() << std::endl;
            }

            return str.str();
        }

        inline Float getAvgKappa() {
            Float res = 0.0f;

            for ( size_t h = 0; h < nComponents(); ++h ) {
                res += c[ h ].m_kappa;
            }

            return (res > 0.0f ? res / (Float) nComponents() : 0.0f );
        }

        inline Float getMaxKappa() {
            Float res = 0.0f;

            for ( size_t h = 0; h < nComponents(); ++h ) {
                if ( c[ h ].m_kappa > res ) {
                    res = c[ h ].m_kappa;
                }
            }

            return res;
        }

        inline Float getMinKappa() {
            Float res = std::numeric_limits<Float>::infinity();

            for ( size_t h = 0; h < nComponents(); ++h ) {
                if ( c[ h ].m_kappa < res ) {
                    res = c[ h ].m_kappa;
                }
            }

            return res;
        }

        virtual bool isKapaLT( Float val ) {
            for ( size_t h = 0; h < nComponents(); h++ ) {
                if ( c[ h ].m_kappa < val ) {
                    return true;
                }
            }

            return false;
        }

        virtual bool isKapaGT( Float val ) {
            for ( size_t h = 0; h < nComponents(); h++ ) {
                if ( c[ h ].m_kappa > val ) {
                    return true;
                }
            }

            return false;
        }    

        virtual void setInfo( const EmInfo & info ) { m_info.merge( info ); }        
        virtual IStack<Particle> * getParticles() { return &m_info.particles; }


        virtual void release() {
            c.dealloc();
            w.dealloc();
            cdf.dealloc();
            m_info.release();
        }

    public:
        //////////////////////////////////////////////////////////////////////////
        // Caching stuff
        //////////////////////////////////////////////////////////////////////////

        static IMPORTANCE_INLINE Float validityRadius(const VMFM& dist, const Config& config, const float pixel2world = 0.f ) {
            Float result = 0;
            const Float logThr = log(config.cache.klDivThres);
            IMPORTANCE_ASSERT(isReal(logThr));
            Float weights = 0.f;

            //the minimum value of PDF (relative to its global maximum) where can be the record projected during interpolation
            for(int i = 0; i < dist.nComponents(); i++) {

#if 1
                const Float clampedKappa = clamp(1 + logThr/dist.c[i].m_kappa, 1e-2f, 0.9999999f);
                const Float maxAngle = acos(clampedKappa);
#else
                const float k = dist.c[i].m_kappa;
                const float asinhArg = (1 + config.cache.klDivThres)*sinh(k);
                const float asinhTerm = Ff::asinh(asinhArg);
                //const float asinh

                IMPORTANCE_ASSERT((asinhArg-sinh(asinhTerm))/asinhArg < 1e-4f);
                const Float maxAngle = acos(5.f/4 - asinhTerm*asinhTerm/(4*k*k));
#endif

                IMPORTANCE_ASSERT(isReal(maxAngle));
                IMPORTANCE_ASSERT(maxAngle > 0.f && maxAngle < PI/2.f);
                const Float actualDist = tan(maxAngle) * std::max(pixel2world*config.cache.minRadius, dist.c[i].m_harmonicDistance);
                weights += dist.w[i];


                const float cosineFactor = std::max(0.1f, dot(dist.c[i].m_mi, dist.m_localFrame.n));
                //const float cosineFactor = 1.f;


                result += dist.w[i]*actualDist/cosineFactor;

                IMPORTANCE_ASSERT(isReal(result) && result > 0.f);
            }
            IMPORTANCE_ASSERT(abs(weights-1.f) < 1e-4f);
            return result;
        }

        template<class TRecord>
        static IMPORTANCE_INLINE Float weight_step1(const Hit& hit, const TRecord& record, const Config& config) {
            const Float normalFactor = ::sqrt(std::max(ZERO, ONE - dot(record.normal, Vector3(hit.normal))));
            const Float distFactor = (Vector3(hit.position)-record.position).size()*record.invValidRadius();
            const Float res = 1.f/(distFactor+normalFactor+1e-6f) - config.cache.invGeometryError;
            IMPORTANCE_ASSERT(isReal(res));
            return res;
        }

		template<class TRecord>
        static IMPORTANCE_INLINE Float weight_step2(const Hit& hit, const TRecord& record, const Vector3 &diffClosest, const float distDenom, const Config& config) {
            const Float normalFactor = ::sqrt(std::max(ZERO, ONE - dot(record.normal, Vector3(hit.normal))));
			const Vector3 diff = Vector3(hit.position)-record.position;
            const Float distFactor = diff.size();
			const Float colinearFactor = (ONE - dot(diffClosest, diff)) * config.cache.colinearFactorWeight; 
            const Float res = 1.f/(distFactor+normalFactor+colinearFactor+1e-6f) - config.cache.invGeometryError;
            IMPORTANCE_ASSERT(isReal(res));
            return res;
        }  

    private:
        size_t getSampledComponentIndex( Float u ) const;

    public:
        IStack<VMF> c;
        IStack<Float> w;
        IStack<Float> cdf;
        bool m_initialized;        
        size_t m_nComponents;                
        /// For online EM - iteration counter (data counter)
        unsigned int k;

        Frame m_localFrame;
    
        /// fitting information
        Importance::EmInfo m_info;
    };
}

#pragma warning(pop)

#endif