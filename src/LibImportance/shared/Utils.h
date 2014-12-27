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

#include "NumberUtils.h"
#include "matrix.h"
#include "Sse.h"
#include "QuadVector.h"
#include "Particle.h"
#include "Stack.h"


#ifdef IMPORTANCE_SINGLE_PRECISION
    #ifndef IMP_PI
        #define IMP_PI         3.14159265358979323846f
    #endif
    #ifndef INV_PI
        #define INV_PI       0.31830988618379067154f
    #endif
    #ifndef INV_TWOPI
        #define INV_TWOPI    0.15915494309189533577f
    #endif
    #ifndef SQRT_TWO
        #define SQRT_TWO     1.41421356237309504880f
    #endif
    #ifndef INV_SQRT_TWO
        #define INV_SQRT_TWO 0.70710678118654752440f
    #endif
#else
    #ifndef IMP_PI
        #define IMP_PI         3.14159265358979323846
    #endif
    #ifndef INV_PI
        #define INV_PI       0.31830988618379067154
    #endif
    #ifndef INV_TWOPI
        #define INV_TWOPI    0.15915494309189533577
    #endif
    #ifndef SQRT_TWO
        #define SQRT_TWO     1.41421356237309504880
    #endif
    #ifndef INV_SQRT_TWO
        #define INV_SQRT_TWO 0.70710678118654752440
    #endif
#endif

namespace Importance {

    struct DiscreteDistribution {
    public:
	    /// Allocate memory for a distribution with the given number of entries
	    explicit inline DiscreteDistribution(size_t nEntries = 0) {
		    reserve(nEntries);
		    clear();
	    }

	    /// Clear all entries
	    inline void clear() {
		    m_cdf.clear();
		    m_cdf.push_back(0.0f);
		    m_normalized = false;
	    }

	    /// Reserve memory for a certain number of entries
	    inline void reserve(size_t nEntries) {
		    m_cdf.reserve(nEntries+1);
	    }

	    /// Append an entry with the specified discrete probability
	    inline void append(Float pdfValue) {
		    m_cdf.push_back(m_cdf[m_cdf.size()-1] + pdfValue);
	    }

	    /// Return the number of entries so far
	    inline size_t size() const {
		    return m_cdf.size()-1;
	    }

	    /// Access an entry by its index
	    inline Float operator[](size_t entry) const {
		    return m_cdf[entry+1] - m_cdf[entry];
	    }

	    /// Have the probability densities been normalized?
	    inline bool isNormalized() const {
		    return m_normalized;
	    }

	    /**
	     * \brief Return the original (unnormalized) sum of all PDF entries
	     *
	     * This assumes that \ref normalize() has previously been called
	     */
	    inline Float getSum() const {
		    return m_sum;
	    }

	    /**
	     * \brief Return the normalization factor (i.e. the inverse of \ref getSum())
	     *
	     * This assumes that \ref normalize() has previously been called
	     */
	    inline Float getNormalization() const {
		    return m_normalization;
	    }

	    /**
	     * \brief Normalize the distribution
	     *
	     * Throws an exception when no entries were previously
	     * added to the distribution.
	     *
	     * \return Sum of the (previously unnormalized) entries
	     */
	    inline Float normalize() {
		    IMPORTANCE_ASSERT(m_cdf.size() > 1);
		    m_sum = m_cdf[m_cdf.size()-1];
		    if (m_sum > 0) {
			    m_normalization = 1.0f / m_sum;
			    for (size_t i=1; i<m_cdf.size(); ++i)
				    m_cdf[i] *= m_normalization;
			    m_cdf[m_cdf.size()-1] = 1.0f;
			    m_normalized = true;
		    } else {
			    m_normalization = 0.0f;
		    }
		    return m_sum;
	    }

	    /**
	     * \brief %Transform a uniformly distributed sample to the stored distribution
	     *
	     * \param[in] sampleValue
	     *     An uniformly distributed sample on [0,1]
	     * \return
	     *     The discrete index associated with the sample
	     */
	    inline size_t sample(Float sampleValue) const {
		    std::vector<Float>::const_iterator entry =
				    std::lower_bound(m_cdf.begin(), m_cdf.end(), sampleValue);
		    size_t index = std::min(m_cdf.size()-2,
			    (size_t) std::max((ptrdiff_t) 0, entry - m_cdf.begin() - 1));

		    /* Handle a rare corner-case where a entry has probability 0
		       but is sampled nonetheless */
		    while (operator[](index) == 0 && index < m_cdf.size()-1)
			    ++index;

		    return index;
	    }

	    /**
	     * \brief %Transform a uniformly distributed sample to the stored distribution
	     *
	     * \param[in] sampleValue
	     *     An uniformly distributed sample on [0,1]
	     * \param[out] pdf
	     *     Probability value of the sample
	     * \return
	     *     The discrete index associated with the sample
	     */
	    inline size_t sample(Float sampleValue, Float &pdf) const {
		    size_t index = sample(sampleValue);
		    pdf = operator[](index);
		    return index;
	    }

	    /**
	     * \brief %Transform a uniformly distributed sample to the stored distribution
	     *
	     * The original sample is value adjusted so that it can be "reused".
	     *
	     * \param[in, out] sampleValue
	     *     An uniformly distributed sample on [0,1]
	     * \return
	     *     The discrete index associated with the sample
	     */
	    inline size_t sampleReuse(Float &sampleValue) const {
		    size_t index = sample(sampleValue);
		    sampleValue = (sampleValue - m_cdf[index])
			    / (m_cdf[index + 1] - m_cdf[index]);
		    return index;
	    }

	    /**
	     * \brief %Transform a uniformly distributed sample.
	     *
	     * The original sample is value adjusted so that it can be "reused".
	     *
	     * \param[in,out]
	     *     An uniformly distributed sample on [0,1]
	     * \param[out] pdf
	     *     Probability value of the sample
	     * \return
	     *     The discrete index associated with the sample
	     */
	    inline size_t sampleReuse(Float &sampleValue, Float &pdf) const {
		    size_t index = sample(sampleValue, pdf);
		    sampleValue = (sampleValue - m_cdf[index])
			    / (m_cdf[index + 1] - m_cdf[index]);
		    return index;
	    }

        inline void release() {
            std::vector<Float>().swap( m_cdf );
        }

	    /**
	     * \brief Turn the underlying distribution into a
	     * human-readable string format
	     */
	    std::string toString() const {
		    std::ostringstream oss;
		    oss << "DiscreteDistribution[sum=" << m_sum << ", normalized="
			    << (int) m_normalized << ", cdf={";
		    for (size_t i=0; i<m_cdf.size(); ++i) {
			    oss << m_cdf[i];
			    if (i != m_cdf.size()-1)
				    oss << ", ";
		    }
		    oss << "}]";
		    return oss.str();
	    }
    private:
	    std::vector<Float> m_cdf;
	    Float m_sum, m_normalization;
	    bool m_normalized;
    };

    /// \brief Tests 2 floats for equivalence with epsilon tolerance
    /// \return true if abs(number1, number2) < 1e-4
    IMPORTANCE_INLINE bool epsEqual(const float number1, const float number2, const float epsilon = 1e-4f) {
        return abs(number1 - number2) < epsilon;
    }

    // udela nahodnej vyber podle random a zapise do random novou hodnotu randomu
    template<class TIter, class TFloat>
    IMPORTANCE_INLINE int discreteSelection(TFloat& random, const TIter cdfBegin, const TIter cdfEnd) {
        const auto searched = (random * cdfEnd[-1]);
        TIter res = std::upper_bound(cdfBegin, cdfEnd, searched);
        if(res == cdfEnd) {
            --res;  //TODO - zalogovat, chybova hlaska?
        }
        const TFloat prevCdf = (res == cdfBegin) ? 0.f : res[-1];
        random = (searched - prevCdf)/(*res - prevCdf);
        IMPORTANCE_ASSERT(isReal(random) && random >= 0.f && random <= 1.f);
        return int(res - cdfBegin);
    }


    IMPORTANCE_INLINE Vector3 squareToHemisphere( Float u1, Float u2 ) {
        Float z = u2;
        Float tmp = std::sqrt(std::max((Float) 0, 1.0f - z*z));        
        Float phi = 2.0f * IMP_PI * u1;
        return Vector3( cos( phi ) * tmp, sin( phi ) * tmp, z );
    }

    // problem with singularity!!!!
    IMPORTANCE_INLINE void hemisphereToSquare( const Vector3 & lDir, Float & u, Float & v ) {       
        IMPORTANCE_ASSERT( lDir.isNormalized() );              
        u = lDir.z;
        Float tmp = std::sqrt(std::max((Float) 0.f, 1.0f - lDir.z*lDir.z));
        Float sinPhi = lDir.y / tmp;
        Float phi = std::asinf( sinPhi );
        v = phi / ( 2.0f * IMP_PI );        
    }

    IMPORTANCE_INLINE Vector3 squareToSphere( Float u1, Float u2 ) {
        Float z = 1.0f - 2.0f * u2;
        Float r = std::sqrt(std::max((Float) 0.0f, 1.0f - z*z));
        Float phi = 2.0f * IMP_PI * u1;        
        return Vector3(r * cos( phi ), r * sin( phi ), z);
    }


    IMPORTANCE_INLINE void coordinateSystem( const Vector3 &a, Vector3 &b, Vector3 &c ) {
        if (std::abs(a.y) < 0.99f) {
            Float invLen = 1.0f / std::sqrt(a.x * a.x + a.z * a.z);
            c = Vector3(a.z * invLen, 0.0f, -a.x * invLen);
        } else {
            Float invLen = 1.0f / std::sqrt(a.y * a.y + a.z * a.z);
            c = Vector3(0.0f, a.z * invLen, -a.y * invLen);
        }
        b = cross(c, a);

        IMPORTANCE_ASSERT(a.isNormalized() && b.isNormalized() && c.isNormalized() && abs(dot(a, b)) < 1e-4f && abs(dot(a, c)) < 1e-4f);
    }

    IMPORTANCE_INLINE Float squareToUniformConePdf ( Float cosCutoff ) {
        return INV_TWOPI / ( 1.f - cosCutoff );
    }

    /// Square root variant that gracefully handles arguments < 0 that are due to roundoff errors
    IMPORTANCE_INLINE float safe_sqrt( float value ) {
        return std::sqrt( std::max( 0.0f, value ) );
    }

    IMPORTANCE_INLINE Float4 safe_sqrt( Float4 value ) {
        return Float4::max(value, Float4(0.f)).sqrt();
    }

    IMPORTANCE_INLINE double safe_sqrt( double value ) {
        return std::sqrt( std::max( 0.0, value ) );
    }

	/// Arccosine variant that gracefully handles arguments > 1 that are due to roundoff errors
	IMPORTANCE_INLINE float safe_acos( float value ) {
		return std::acos( std::min( 1.0f, std::max( -1.0f, value ) ) );
	}

	/// Arccosine variant that gracefully handles arguments > 1 that are due to roundoff errors
	IMPORTANCE_INLINE double safe_acos( double value ) {
		return std::acos( std::min( 1.0, std::max( -1.0, value ) ) );
	}

    //IMPORTANCE_INLINE void sincos( float theta, float *_sin, float *_cos ) {
    //    *_sin = sinf(theta);
    //    *_cos = cosf(theta);
    //}

    //IMPORTANCE_INLINE void sincos( double theta, double *_sin, double *_cos ) {
    //    *_sin = sin(theta);
    //    *_cos = cos(theta);
    //}

    IMPORTANCE_INLINE Vector3 squareToUniformCone( Float cosCutoff, const Vector2 & sample ) {
        Float cosTheta = ( 1.0f - sample.x ) + sample.x * cosCutoff;
        Float sinTheta = safe_sqrt( 1.0f - cosTheta * cosTheta );

        Float sinPhi, cosPhi;
        Ff::sincos( 2.0f * IMP_PI * sample.y, sinPhi, cosPhi );

        return Vector3(cosPhi * sinTheta,
            sinPhi * sinTheta, cosTheta);
    }

    IMPORTANCE_INLINE Vector2 squareToUniformDiskConcentric( const Vector2 & sample ) {
        Float r1 = 2.0f*sample.x - 1.0f;
        Float r2 = 2.0f*sample.y - 1.0f;

        Vector2 coords;
        if ( r1 == 0 && r2 == 0 ) {
            coords = Vector2( 0, 0 );
        } else if ( r1 > -r2 ) { /* Regions 1/2 */
            if ( r1 > r2 )
                coords = Vector2( r1, (IMP_PI/4.0f) * r2/r1 );
            else
                coords = Vector2( r2, (IMP_PI/4.0f) * (2.0f - r1/r2) );
        } else { /* Regions 3/4 */
            if ( r1 < r2 )
                coords = Vector2( -r1, (IMP_PI/4.0f) * (4.0f + r2/r1) );
            else
                coords = Vector2( -r2, (IMP_PI/4.0f) * (6.0f - r1/r2) );
        }

        Vector2 result;
        Ff::sincos( coords.y, result.y, result.x );
        result.x *= coords.x; 
        result.y *= coords.x;
        IMPORTANCE_ASSERT(isReal(result.x) && isReal(result.y));
        return result;
    }

    IMPORTANCE_INLINE Vector3 mapDirection(const int x, const int y, const int width, const int height) {
        const float phi = (2*PI*x)/width;
        const float theta = float(y)/height;
        const float cosTheta = 2*theta - 1;
        const float sinTheta = sqrt(1-cosTheta*cosTheta);
        float cosPhi, sinPhi;
        Ff::sincos(phi, sinPhi, cosPhi);
        const Vector3 result(cosPhi*sinTheta, sinPhi*sinTheta, cosTheta);
        IMPORTANCE_ASSERT(result.isNormalized());
        return result;
    }



    /// \brief Transforms sample from unit square to uniformly sampled hemisphere using shirley's mapping
    IMPORTANCE_INLINE Vector3 uniformHemisphereShirley( Float u1, Float u2 ) {
        Float phi, r;
        const Float a = 2.f*u1 - 1.f;   // (a,b) is now on [-1,1]^2
        const Float b = 2.f*u2 - 1.f;
        if(a > -b) {     // region 1 or 2 
            if(a > b) {  // region 1, also |a| > |b|
                r = a;
                phi = (PI/4) * (b/a);
            } else {  // region 2, also |b| > |a|
                r = b;
                phi = (PI/4) * (2.f - (a/b));
            }
        } else {        // region 3 or 4 
            if (a < b) {  // region 3, also |a| >= |b|, a != 0 
                r = -a;
                phi = (PI/4) * (4 + (b/a));
            }
            else       {  // region 4, |b| >= |a|, but a==0 and b==0 could occur.
                r = -b;
                if (b != 0.f) {
                    phi = (PI/4) * (6.f - (a/b));
                } else {
                    phi = 0.f;
                }
            }
        }
        const Float rSqr = r*r;
        const Float scaleFactor = r * safe_sqrt( 2 - rSqr );
        IMPORTANCE_ASSERT( isReal( r + phi ) );
        Float cosPhi, sinPhi;
        Ff::sincos( phi, sinPhi, cosPhi );
        const Vector3 result( cosPhi * scaleFactor, sinPhi * scaleFactor, 1 - rSqr );
        IMPORTANCE_ASSERT( result.isNormalized() );
        return result;
    }

    IMPORTANCE_INLINE QuadVector3 uniformHemisphereShirley( Float4 u1, Float4 u2 ) {
        //Float4 phi, r;
        const Float4 a = 2.f*u1 - 1.f;   // (a,b) is now on [-1,1]^2
        const Float4 b = 2.f*u2 - 1.f;

        const Float4 r = (a > -b).blend(Float4::max(a, b), -Float4::min(a, b));

        const Float4 aOverB = a/b;
        const Float4 bOverA = aOverB.getInverse();
        Float4 phi = (a > -b).blend(
                (a > b).blend(bOverA, 2-aOverB),
                (a < b).blend(4+bOverA, (b != Float4(0.f)).blend(6-aOverB, Float4(0.f)))
            );
        phi *= PI/4;

        const Float4 rSqr = r*r;
        const Float4 scaleFactor = r * safe_sqrt( 2 - rSqr );
        IMPORTANCE_ASSERT( isReal( r + phi ) );
        Float4 cosPhi, sinPhi;
        Ff::sincos( phi.data, sinPhi.data, cosPhi.data );
        const QuadVector3 result( cosPhi * scaleFactor, sinPhi * scaleFactor, 1 - rSqr );
        IMPORTANCE_ASSERT( result.isNormalized() );
#ifdef LIBIMP_DEBUG
        for(int i = 0; i < 4; ++i) {
            const Vector3 thisResult = result.get(i);
            const Vector3 otherResult = uniformHemisphereShirley(u1[i], u2[i]);
            IMPORTANCE_ASSERT((thisResult-otherResult).size() < 1e-3f);
        }
#endif
        return result;
    }

    /// \brief The mapping from unit z-up hemisphere to unit square
    IMPORTANCE_INLINE void uniformHemisphereShirleyInverse( const Vector3 dir, Float & _u1, Float & _u2 ) {
        const Float r = safe_sqrt( 1.f - dir.z );
        Float phi = std::atan2( dir.y, dir.x );

        Float a, b;
        if (phi < -PI/4) {
            phi += 2*PI;        // in range [-pi/4,7pi/4] 
        }
        if ( phi < PI/4) {                   // region 1 
            a = r;
            b = phi * a / (PI/4);
        } else if ( phi < 3*PI/4 ) {         // region 2 
            b = r;
            a = -(phi - PI/2) * b / (PI/4);
        } else if ( phi < 5*PI/4 ) {         // region 3
            a = -r;
            b = (phi - PI) * a / (PI/4);
        } else {                             // region 4
            b = -r;
            a = -(phi - 3*PI/2) * b / (PI/4);
        }

        _u1 = ( a + 1.f ) * 0.5f;
        _u2 = ( b + 1.f ) * 0.5f;
#ifdef LIBIMP_DEBUG
        const Vector3 temp = uniformHemisphereShirley(_u1, _u2);
        IMPORTANCE_ASSERT(abs(dir.x - temp.x) < 1e-3f);
        IMPORTANCE_ASSERT(abs(dir.y - temp.y) < 1e-3f);
        IMPORTANCE_ASSERT(abs(dir.z - temp.z) < 1e-3f);
#endif
    }

    IMPORTANCE_INLINE void uniformHemisphereShirleyInverse( const QuadVector3& dir, Float4 & _u1, Float4 & _u2 ) {
        const Float4 r = safe_sqrt( Float4(1.f) - dir.z );
        Float4 phi = Float4(Ff::atan2( dir.y.data, dir.x.data));

        phi = (phi < -PI/4).blend(phi+2*PI, phi);

        Float4 a, b;

        for(int i = 0; i < 4; ++i) {
            if ( phi[i] < PI/4) {                   // region 1 
                a[i] = r[i];
                b[i] = phi[i] * a[i] / (PI/4);
            } else if ( phi[i] < 3*PI/4 ) {         // region 2 
                b[i] = r[i];
                a[i] = -(phi[i] - PI/2) * b[i] / (PI/4);
            } else if ( phi[i] < 5*PI/4 ) {         // region 3
                a[i] = -r[i];
                b[i] = (phi[i] - PI) * a[i] / (PI/4);
            } else {                             // region 4
                b[i] = -r[i];
                a[i] = -(phi[i] - 3*PI/2) * b[i] / (PI/4);
            }
        }
        
        _u1 = ( a + Float4(1.f) ) * 0.5f;
        _u2 = ( b + Float4(1.f) ) * 0.5f;
#ifdef LIBIMP_DEBUG
        for(int i = 0; i < 4; ++i) {
            const Vector2 thisResult(_u1[i], _u2[i]);
            Vector2 otherResult;
            uniformHemisphereShirleyInverse(dir.get(i), otherResult.x, otherResult.y);
            IMPORTANCE_ASSERT((thisResult-otherResult).size() < 1e-3f);
        }
#endif
    }


    IMPORTANCE_INLINE Vector3 uniformSphere(const float rot, const float elev) {
        const float cosTheta = 2*elev - 1;
        const float sinTheta = sqrt(1-cosTheta*cosTheta);
        const float phi = 2*PI*rot;
        float cosPhi, sinPhi;
        Ff::sincos(phi, sinPhi, cosPhi);
        const Vector3 result(cosPhi*sinTheta, sinPhi*sinTheta, cosTheta);
        IMPORTANCE_ASSERT(result.isNormalized());
        return result;
    }
    IMPORTANCE_INLINE void uniformSphereInverse(const Vector3 dir, float& rot, float& elev) {
        IMPORTANCE_ASSERT(dir.isNormalized());
        const float cosTheta = dir.z;
        Vector2 result;
        elev = (cosTheta+1)/2;
        float phi = Ff::atan2(dir.y, dir.x);
        if(phi < 0.f) {
            phi += 2*PI;
        }
        rot = phi / (2*PI);
#ifdef LIBIMP_DEBUG
        const Vector3 check = uniformSphere( rot, elev );
        IMPORTANCE_ASSERT(epsEqual(check.x, dir.x));
        IMPORTANCE_ASSERT(epsEqual(check.y, dir.y));
        IMPORTANCE_ASSERT(epsEqual(check.z, dir.z));
#endif
    }

    IMPORTANCE_INLINE Vector2 hemisphereToDisk(const Vector3 in) {
        IMPORTANCE_ASSERT( in.isNormalized() );
        Vector2 result(in.x, in.y);
        const float radius = std::sqrt(1-in.z);        
        if ( radius > 0.f ) {
            result /= result.size();
        }
        result *= radius;
        //IMPORTANCE_ASSERT(std::min(result.x, result.y) >= -1.f && std::max(result.x, result.y) <= 1.f);
        IMPORTANCE_ASSERT( result.isValid() );
        return result;
    }
    IMPORTANCE_INLINE Vector3 diskToHemisphere(const Vector2 in) {        
        const float z = 1.f - in.square();
        IMPORTANCE_ASSERT( z >= 0.f );
        /* Distance of projection on the floor. */
        const float radius = std::sqrt(1.f-z*z);
        /* Carefully normalize the direction on the floor. */
        Vector2 dir( in );        
        if( radius > 0.f ) {
            dir /= in.size();
        }
        dir *= radius;
        const Vector3 result(dir.x, dir.y, z);
        IMPORTANCE_ASSERT(result.isNormalized());
        
//#ifdef LIBIMP_DEBUG
#if 0
        const Vector2 temp = hemisphereToDisk(result);
        IMPORTANCE_ASSERT(abs(temp.x-in.x) < 1e-3f);
        IMPORTANCE_ASSERT(abs(temp.y-in.y) < 1e-3f
            );
#endif
        return result;
    }

    inline Vector3 squareToCosineHemisphere( const Vector2 &sample ) {
        Vector2 p = squareToUniformDiskConcentric( sample );
        Float z = std::sqrt( std::max( 0.f, 1.0f - p.x*p.x - p.y*p.y ) );

        /* Guard against numerical imprecisions */
        if ( z == 0 )
            z = 1e-10f;

        return Vector3( p.x, p.y, z );
    }

    inline Float squareToCosineHemispherePdf( const Vector3 & d ) { 
        // d must be in a local coordinate system !!        
        return INV_PI * d.z;
    }

    template<class TIterator>
    IMPORTANCE_INLINE Float avgDistance( TIterator begin, TIterator end ) {
        //IMPORTANCE_ASSERT( begin <= end );        
        size_t n = end - begin;

        if ( n == 0 ) {
            return 1.f;
        }

        Float res = 0.f;
        while ( begin != end ) {
            res += begin->getDistance();
            ++begin;
        }

        IMPORTANCE_ASSERT( res > 0.f );
        return n / res;
    }


    template<class TIterator>
    IMPORTANCE_INLINE Float4 avgDistanceV( TIterator begin, TIterator end ) {
        IMPORTANCE_ASSERT( begin <= end );        
        size_t n = end - begin;

        if ( n == 0 ) {
            return Float4(1.f);
        }

        Float4 res = Float4(0.f);
        while ( begin != end ) {
            res += begin->getDistance();
            ++begin;
        }

        return Float4((float)n) / res;
    }

    template<class TIterator>
    IMPORTANCE_INLINE Float harmonicDistanceFromInverseDistances( TIterator begin, TIterator end ) {
        //IMPORTANCE_ASSERT( begin <= end );        
        size_t n = end - begin;

        if ( n == 0 ) {
            return 1.f;
        }

        Float res = 0.f;
        while ( begin != end ) {
            res += begin->getInverseDistance();
            ++begin;
        }

        IMPORTANCE_ASSERT( res > 0.f );
        return n / res;
    }

    inline float linePointDistance(const Vector3 lineO, const Vector3 lineD, const Vector3 point) {
        return ((lineO-point)-(dot(lineO-point, lineD)*lineD)).size();
    }

    namespace Transform {
        IMPORTANCE_INLINE Matrix3x3 rotation(const Vector3 & naxis, Float angle) {
            Float sinTheta, cosTheta;

            /* Make sure that the axis is normalized */
            IMPORTANCE_ASSERT( naxis.isNormalized() );
            Ff::sincos(angle, sinTheta, cosTheta);

            Matrix3x3 result;
            result(0, 0) = naxis.x * naxis.x + (1.0f - naxis.x * naxis.x) * cosTheta;
            result(0, 1) = naxis.x * naxis.y * (1.0f - cosTheta) - naxis.z * sinTheta;
            result(0, 2) = naxis.x * naxis.z * (1.0f - cosTheta) + naxis.y * sinTheta;

            result(1, 0) = naxis.x * naxis.y * (1.0f - cosTheta) + naxis.z * sinTheta;
            result(1, 1) = naxis.y * naxis.y + (1.0f - naxis.y * naxis.y) * cosTheta;
            result(1, 2) = naxis.y * naxis.z * (1.0f - cosTheta) - naxis.x * sinTheta;

            result(2, 0) = naxis.x * naxis.z * (1.0f - cosTheta) - naxis.y * sinTheta;
            result(2, 1) = naxis.y * naxis.z * (1.0f - cosTheta) + naxis.x * sinTheta;
            result(2, 2) = naxis.z * naxis.z + (1.0f - naxis.z * naxis.z) * cosTheta;
            
            return result;
        }
   }   

    template<class T>
    IMPORTANCE_INLINE T lerp(const T& a, const T& b, const float amountB) {
        return a*(1-amountB) + b*amountB;
    }


    struct SampleStats {
        double var;
        double mean;
        size_t count;
        double max, min;
        //Float modus; // computationally too expensive - use 3rd part library if you really want this
        double median;

        SampleStats() {
            var = mean = max = min = median = -1.f;
            count = -1;
        }
    
        std::string toString() const {
            std::ostringstream ostr;
            ostr << "[" << std::endl
                 << "mean: " << mean << std::endl
                 << "var: " << var <<std::endl
                 << "median: " << median << std::endl
                 << "count: " << count << std::endl
                 << "max: " << max << std::endl
                 << "min: " << min << std::endl
                 << "]" << std::endl;
            return ostr.str();
        }
    };

    class ParticlesComparator {
    private:
        const IStack<Particle> & m_p;
    public:
        ParticlesComparator( const IStack<Particle> & particles ) : m_p( particles ) {}
        bool operator()( const int & a, const int & b ) const {
            IMPORTANCE_ASSERT( a >= 0 && b >= 0 && a < m_p.size() && b < m_p.size() );
            return m_p[ a ].weight < m_p[ b ].weight;
        }
        ParticlesComparator( const ParticlesComparator & pc ) : m_p( pc.m_p ) {}
        ParticlesComparator & operator=( const ParticlesComparator & pc ) { 
            this->ParticlesComparator::ParticlesComparator( pc );
        }
    };

    IMPORTANCE_INLINE void computeSampleStatistics( const IStack<Particle> & samples, SampleStats & stats ) {
        stats.mean = 0.0;
        stats.count = samples.size();
        stats.max = -std::numeric_limits<double>::max();
        stats.min = std::numeric_limits<double>::max();

        if ( stats.count == 0 ) {
            return;
        }

        IStack<int> indices( stats.count );
        int i = 0;
        for ( auto it = samples.cbegin(); it < samples.cend(); ++it, ++i ) {
            stats.mean += it->weight;
            stats.max = std::max( stats.max, (double) it->weight );
            stats.min = std::min( stats.min, (double) it->weight );
            indices[ i ] = i;
        }
        stats.mean /= stats.count;

        stats.var = 0.f;
        for ( auto it = samples.cbegin(); it < samples.cend(); ++it ) {
            stats.var += ( ( stats.mean - it->weight ) * ( stats.mean - it->weight ) );
        }
        stats.var /= (double) std::max( (size_t) 1, stats.count - 1 );

        std::nth_element( indices.begin(), indices.begin() + indices.size() / 2, indices.end(),
            ParticlesComparator( samples ) );
        stats.median = samples[ indices[ indices.size() / 2 ] ].weight; 
    }
}

#include "frame.h"

namespace Importance {
    class SquareHemisphereMapping {
    public:
        IMPORTANCE_INLINE bool fromSquare( const Frame & coords, const Vector2 & point, Vector3 & oWorldDir ) const {
            if ( point.x < 0.0f || point.x > 1.0f || point.y < 0.0f || point.y > 1.0f ) {
                oWorldDir = coords.toWorld( Vector3(0.0f,0.0f,-1.0f) );
                return false;
            }
            const Vector3 localDir = uniformHemisphereShirley( point.x, point.y );
            oWorldDir              =  coords.toWorld( localDir );
            IMPORTANCE_ASSERT( oWorldDir.isReal() );
            return true;
        }

        IMPORTANCE_INLINE bool fromSquare( const Vector2 & point, Vector3 & oLocalDir ) const {
            if ( point.x < 0.0f || point.x > 1.0f || point.y < 0.0f || point.y > 1.0f ) {
                oLocalDir = Vector3(0.0f,0.0f,-1.0f);
                return false;
            }
            oLocalDir = uniformHemisphereShirley( point.x, point.y );            
            IMPORTANCE_ASSERT( oLocalDir.isReal() );
            return true;
        }

        IMPORTANCE_INLINE bool toSquare( const Frame & coords, const Vector3 & worldDir, Vector2 & oPoint ) const {
            const Vector3 localDir = coords.toLocal( worldDir );
            if ( localDir.z < 0.0f ) {
                return false;
            }
            uniformHemisphereShirleyInverse( localDir, oPoint.x, oPoint.y );
            return true;
        }

        IMPORTANCE_INLINE bool toSquare( const Vector3 & localDir, Vector2 & oPoint ) const {
            if ( localDir.z < 0.0f ) {
                return false;
            }
            uniformHemisphereShirleyInverse( localDir, oPoint.x, oPoint.y );
            return true;
        }

        static IMPORTANCE_INLINE Float jacobian() {
            return 2.0f * IMP_PI;
        }
    };


    class SquareSphereMapping {
    public:
        IMPORTANCE_INLINE bool fromSquare( const Frame & /*coords*/, const Vector2 & point, Vector3 & oWorldDir ) const {
            return fromSquare( point, oWorldDir );
        }

        IMPORTANCE_INLINE bool fromSquare( const Vector2 & point, Vector3 & oWorldDir ) const {
            oWorldDir = uniformSphere( point.x, point.y );
            return true;
        }

        IMPORTANCE_INLINE bool toSquare( const Frame & /*coords*/, const Vector3 & worldDir, Vector2 & oPoint ) const {
            return toSquare( worldDir, oPoint );
        }

        IMPORTANCE_INLINE bool toSquare( const Vector3 & worldDir, Vector2 & oPoint ) const {
            uniformSphereInverse( worldDir, oPoint.x, oPoint.y );
            return true;
        }

        static IMPORTANCE_INLINE Float jacobian() {
            return 4.0f * IMP_PI;
        }
    };

}