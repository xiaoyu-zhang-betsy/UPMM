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
#include "../shared/FastFloat/FastFloat.h"
#include "../shared/StaticStack.h"
#include "../shared/utils.h"
#include "../shared/QuadVector.h"
#include "../shared/matrix.h"
#include "../shared/Vector2.h"
#include "../vmf/random.h"
#include "../shared/frame.h"
#include "../caching/CacheStats.h"
#include "arraysse.h"
#include "../distribution/DefaultDistributionModel.h"

#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>


#pragma warning(push)
#pragma warning(disable:4265)
#pragma warning(disable:4127)
#pragma warning(disable:4100)


namespace Importance {

    typedef Vector2t<VectorTraits< Float4 > > QuadVector2;

	static const double a[] =
	{
		-3.969683028665376e+01,
		2.209460984245205e+02,
		-2.759285104469687e+02,
		1.383577518672690e+02,
		-3.066479806614716e+01,
		2.506628277459239e+00
	};

	static const double b[] =
	{
		-5.447609879822406e+01,
		1.615858368580409e+02,
		-1.556989798598866e+02,
		6.680131188771972e+01,
		-1.328068155288572e+01
	};

	static const double c[] =
	{
		-7.784894002430293e-03,
		-3.223964580411365e-01,
		-2.400758277161838e+00,
		-2.549732539343734e+00,
		4.374664141464968e+00,
		2.938163982698783e+00
	};

	static const double d[] =
	{
		7.784695709041462e-03,
		3.224671290700398e-01,
		2.445134137142996e+00,
		3.754408661907416e+00
	};

#define LOW 0.02425
#define HIGH 0.97575
	inline double inverseGaussianCDF(double p)
	{
		double q, r;

		errno = 0;

		if (p < 0 || p > 1)
		{
			errno = EDOM;
			return 0.0;
		}
		else if (p == 0)
		{
			errno = ERANGE;
			return -HUGE_VAL /* minus "infinity" */;
		}
		else if (p == 1)
		{
			errno = ERANGE;
			return HUGE_VAL /* "infinity" */;
		}
		else if (p < LOW)
		{
			/* Rational approximation for lower region */
			q = sqrt(-2 * log(p));
			return (((((c[0] * q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
				((((d[0] * q + d[1])*q + d[2])*q + d[3])*q + 1);
		}
		else if (p > HIGH)
		{
			/* Rational approximation for upper region */
			q = sqrt(-2 * log(1 - p));
			return -(((((c[0] * q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
				((((d[0] * q + d[1])*q + d[2])*q + d[3])*q + 1);
		}
		else
		{
			/* Rational approximation for central region */
			q = p - 0.5;
			r = q*q;
			return (((((a[0] * r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
				(((((b[0] * r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1);
		}
	}


	const double gaussianConstant[6] = { 0.2316419, 0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429 };
	inline double normalDistribution(double x){
		return exp(-0.5 * x * x) / sqrt(2.0 * IMP_PI);
	}
	inline double gaussianCDF(double x){
		bool negx = false;
		if (x < 0.0){
			x = -x;
			negx = true;
		}
		double phi = normalDistribution(x);
		double t = 1.0 / (1.0 + gaussianConstant[0] * x);
		double tsum = 0.0;
		double ti = 1.0;
		for (int i = 1; i <= 5; i++){
			ti *= t;
			tsum += ti * gaussianConstant[i];
		}
		double res = (negx) ? phi * tsum : 1.0 - phi * tsum;
		res = std::min(1.0, std::max(0.0, res));
		return res;
	}	

    template<typename TScalar>
    class GaussianSamplingLobes {
    public:
        typedef TScalar ScalarType;
        static const int WIDTH = 4;

        QuadVector2 mean;
        
        // precision matrix!!! (i.e. inverse of covariance)
        SPDMatrix2x2< TScalar > cov;
        
        TScalar normalization;
        TScalar weights;
    
        IMPORTANCE_INLINE TScalar pdf(const QuadVector2 & v) const {
            const TScalar e = TScalar( Ff::importanceExp( (TScalar(-0.5f) * cov.quadraticForm( v - mean ) ).data ) );
            const TScalar res = e * normalization;
            IMPORTANCE_ASSERT( res.isReal() );
            return res;
        }

        IMPORTANCE_INLINE TScalar pdf(const TScalar & x, const TScalar & y) const {
            return pdf( QuadVector2(x,y) );
        }

        IMPORTANCE_INLINE Vector2 sampleSquare(const int index, Vector2 random) const {
            IMPORTANCE_ASSERT( random.x >= 0.f && random.x <= 1.f );
            // First sample x from two normal distributions using Box-Muller method
            random.x = std::max( random.x, std::numeric_limits<Float>::min() );
            const Float mult = std::sqrt( - 2.0f * Ff::log( random.x ) );
            const Float angle = 2 * IMP_PI * random.y; // TODO: test -pi to pi sampling
            Float sin_angle, cos_angle;
            Ff::sincos( angle, sin_angle, cos_angle );
            const Vector2 x = Vector2( mult * cos_angle, mult * sin_angle );
            // Once I have vector x from N(0,1) I can transform it to v = A * x + mean, where A is matrix from equation Covariance = A * A^TScalar
            float c2 = cov.m[ 1 ][ index ] * cov.m[ 1 ][ index ];
            float invDetSqrt = 1.0f / std::sqrt( cov.m[ 0 ][ index ] * cov.m[ 2 ][ index ] - c2 );
            Vector2 dir;
            if ( cov.m[0][ index ] > cov.m[2][ index ] )
            {
                float a22 = std::sqrt( cov.m[ 0 ][ index ] );
                float a12 = - cov.m[ 1 ][ index ] / a22;
                float a11 = std::sqrt( cov.m[ 2 ][ index ] - c2 / cov.m[ 0 ][ index ] );
                dir.x = invDetSqrt * ( a11 * x.x + a12 * x.y ) + mean.x[ index ];
                dir.y = invDetSqrt * a22 * x.y + mean.y[ index ];
            }
            else
            {
                float a11 = std::sqrt( cov.m[ 2 ][ index ] );
                float a21 = - cov.m[ 1 ][ index ] / a11;
                float a22 = std::sqrt( cov.m[ 0 ][ index ] - c2 / cov.m[ 2 ][ index ] );
                dir.x = invDetSqrt * a11 * x.x + mean.x[ index ];
                dir.y = invDetSqrt * ( a21 * x.x + a22 * x.y ) + mean.y[ index ];
            }            

            return dir;
        }

        IMPORTANCE_INLINE Vector2 sample(const int index, Vector2 random) const {
            return sampleSquare( index, random );            
        }

        IMPORTANCE_INLINE void lightweightInit() {
            const TScalar det_cov = cov.determinant();
            IMPORTANCE_ASSERT(det_cov.isReal());
            normalization = det_cov.abs().sqrt() * 1/ (  2.0f * IMP_PI);
            IMPORTANCE_ASSERT( normalization.isReal() );
            IMPORTANCE_ASSERT( cov.m[ 0 ].isReal() );
            IMPORTANCE_ASSERT( cov.m[ 1 ].isReal() );
            IMPORTANCE_ASSERT( cov.m[ 2 ].isReal() );
        }

		IMPORTANCE_INLINE Vector2 toLobe(int index, Vector2 dir) const{
			float c2 = cov.m[1][index] * cov.m[1][index];
			float invDetSqrt = 1.0f / std::sqrt(cov.m[0][index] * cov.m[2][index] - c2);
			Vector2 x;
			if (cov.m[0][index] > cov.m[2][index])
			{
				float a22 = std::sqrt(cov.m[0][index]);
				float a12 = -cov.m[1][index] / a22;
				float a11 = std::sqrt(cov.m[2][index] - c2 / cov.m[0][index]);
				// 				dir.x = invDetSqrt * (a11 * x.x + a12 * x.y) + mean.x[index];
				// 				dir.y = invDetSqrt * a22 * x.y + mean.y[index];
				x.y = (dir.y - mean.y[index]) / (invDetSqrt * a22);
				x.x = ((dir.x - mean.x[index]) / invDetSqrt - a12 * x.y) / a11;
			}
			else
			{
				float a11 = std::sqrt(cov.m[2][index]);
				float a21 = -cov.m[1][index] / a11;
				float a22 = std::sqrt(cov.m[0][index] - c2 / cov.m[2][index]);
// 				dir.x = invDetSqrt * a11 * x.x + mean.x[index];
// 				dir.y = invDetSqrt * (a21 * x.x + a22 * x.y) + mean.y[index];
				x.x = (dir.x - mean.x[index]) / (invDetSqrt * a11);
				x.y = ((dir.y - mean.y[index]) / invDetSqrt - a21 * x.x) / a22;
			}
			return x;
		}

		IMPORTANCE_INLINE Float gatherAreaPdfLobe(int index, Vector2* criticalPoints, int numCriticalPoints, Vector2* componentBounds, int &topComponentBounds) const{
			// uniform sampling, add default bound and return
			if (numCriticalPoints == 0){
				componentBounds[topComponentBounds] = Vector2(0.f, 1.f);
				componentBounds[topComponentBounds + 1] = Vector2(0.f, 1.f);
				topComponentBounds += 2;
				return 1.f;
			}

			// project corner points to inside lobe coords
			Vector2 xmax = Vector2(-10000.f, -10000.f);
			Vector2 xmin = Vector2(10000.f, 10000.f);
// 			for (int i = 0; i < numCriticalPoints; i++){
// 				Vector2 dir = criticalPoints[i];
// 				Vector2 x = toLobe(index, dir);
// 				xmin.x = std::min(xmin.x, x.x);
// 				xmin.y = std::min(xmin.y, x.y);
// 				xmax.x = std::max(xmax.x, x.x);
// 				xmax.y = std::max(xmax.y, x.y);
// 			}
			float c2 = cov.m[1][index] * cov.m[1][index];
			float invDetSqrt = 1.0f / std::sqrt(cov.m[0][index] * cov.m[2][index] - c2);			
			if (cov.m[0][index] > cov.m[2][index])
			{
				float a22 = std::sqrt(cov.m[0][index]);
				float a12 = -cov.m[1][index] / a22;
				float a11 = std::sqrt(cov.m[2][index] - c2 / cov.m[0][index]);
				float f0 = 1.f / (invDetSqrt * a22);
				float f1 = 1.f / invDetSqrt;
				float f2 = 1.f / a11;

				for (int i = 0; i < numCriticalPoints; i++){
					Vector2 x;
					Vector2 dir = criticalPoints[i];
					x.y = (dir.y - mean.y[index]) * f0;
					x.x = ((dir.x - mean.x[index]) * f1 - a12 * x.y) * f2;

					xmin.x = std::min(xmin.x, x.x);
					xmin.y = std::min(xmin.y, x.y);
					xmax.x = std::max(xmax.x, x.x);
					xmax.y = std::max(xmax.y, x.y);
				}
			}
			else
			{
				float a11 = std::sqrt(cov.m[2][index]);
				float a21 = -cov.m[1][index] / a11;
				float a22 = std::sqrt(cov.m[0][index] - c2 / cov.m[2][index]);

				float f0 = 1.f / (invDetSqrt * a11);
				float f1 = 1.f / invDetSqrt;
				float f2 = 1.f / a22;
				for (int i = 0; i < numCriticalPoints; i++){
					Vector2 x;
					Vector2 dir = criticalPoints[i];
					x.x = (dir.x - mean.x[index]) * f0;
					x.y = ((dir.y - mean.y[index]) * f1 - a21 * x.x) *f2;

					xmin.x = std::min(xmin.x, x.x);
					xmin.y = std::min(xmin.y, x.y);
					xmax.x = std::max(xmax.x, x.x);
					xmax.y = std::max(xmax.y, x.y);
				}
			}

			// calculate the cdfs of this bbox and the probability integral
			Float cdfx0 = gaussianCDF(xmin.x);
			Float cdfx1 = gaussianCDF(xmax.x);
			Float cdfy0 = gaussianCDF(xmin.y);
			Float cdfy1 = gaussianCDF(xmax.y);
			Float prob = (cdfx1 - cdfx0) * (cdfy1 - cdfy0);			
			if (prob <= 0.f){
				componentBounds[topComponentBounds] = Vector2(0.f, 1.f);
				componentBounds[topComponentBounds + 1] = Vector2(0.f, 1.f);
				topComponentBounds += 2;
				return 0.f;
			}
			componentBounds[topComponentBounds] = Vector2(cdfx0, cdfx1);
			componentBounds[topComponentBounds + 1] = Vector2(cdfy0, cdfy1);
			topComponentBounds += 2;
			return prob;
		}

		IMPORTANCE_INLINE Vector2 sampleGatherArea(int index, Vector2 random, 
			Vector2 bound0, Vector2 bound1) const{			
			// First sample x from two normal distributions using Box-Muller method
			random.x = random.x * (bound0.y - bound0.x) + bound0.x;
			random.y = random.y * (bound1.y - bound1.x) + bound1.x;
 			IMPORTANCE_ASSERT(random.x > 0.f && random.x < 1.f);
 			IMPORTANCE_ASSERT(random.y > 0.f && random.y < 1.f);
			Float gx = inverseGaussianCDF(random.x);
			Float gy = inverseGaussianCDF(random.y);
			const Vector2 x = Vector2(gx, gy);			
			// Once I have vector x from N(0,1) I can transform it to v = A * x + mean, where A is matrix from equation Covariance = A * A^TScalar
			float c2 = cov.m[1][index] * cov.m[1][index];
			float invDetSqrt = 1.0f / std::sqrt(cov.m[0][index] * cov.m[2][index] - c2);
			Vector2 dir;
			if (cov.m[0][index] > cov.m[2][index])
			{
				float a22 = std::sqrt(cov.m[0][index]);
				float a12 = -cov.m[1][index] / a22;
				float a11 = std::sqrt(cov.m[2][index] - c2 / cov.m[0][index]);
				dir.x = invDetSqrt * (a11 * x.x + a12 * x.y) + mean.x[index];
				dir.y = invDetSqrt * a22 * x.y + mean.y[index];
			}
			else
			{
				float a11 = std::sqrt(cov.m[2][index]);
				float a21 = -cov.m[1][index] / a11;
				float a22 = std::sqrt(cov.m[0][index] - c2 / cov.m[2][index]);
				dir.x = invDetSqrt * a11 * x.x + mean.x[index];
				dir.y = invDetSqrt * (a21 * x.x + a22 * x.y) + mean.y[index];
			}

			return dir;
		}		
    };

    template<int TLobesCount, typename TLobeType, typename TMapping>
    class GaussianSamplingMix : public DefaultDistributionModel {
    public:
        typedef typename TLobeType::ScalarType TScalar;

        static const int WIDTH = TLobeType::WIDTH;
        ImStaticArray<TLobeType, TLobesCount/TLobeType::WIDTH> lobes;
        Frame localFrame;
        int storedLobes;

		Vector2 fromPolarToSquareThetaPhi(Float costheta, Float sintheta, Float cosphi, Float sinphi) const{
			Float z = costheta;
			Float y = sinphi * sintheta;
			Float x = cosphi * sintheta;
			Vector3 dir(x, y, z);
			Vector2 ret;
			getMapping().toSquare(dir, ret);
			return ret;
		}
		Vector2 fromPolarToSquare(Float theta, Float phi) const{
			Float z = cos(theta);
			Float sintheta = sin(theta);
			Float y = sin(phi) * sintheta;
			Float x = cos(phi) * sintheta;
			Vector3 dir(x, y, z);
			Vector2 ret;
			getMapping().toSquare(dir, ret);
			return ret;
		}
		Vector2 fromPolarToSquareTheta(Float costheta, Float sintheta, Float phi) const{
			Float z = costheta;
			Float y = sin(phi) * sintheta;
			Float x = cos(phi) * sintheta;
			Vector3 dir(x, y, z);
			Vector2 ret;
			getMapping().toSquare(dir, ret);
			return ret;
		}
		Vector2 fromPolarToSquarePhi(Float theta, Float cosphi, Float sinphi) const{
			Float z = cos(theta);
			Float sintheta = sin(theta);
			Float y = sinphi * sintheta;
			Float x = cosphi * sintheta;
			Vector3 dir(x, y, z);
			Vector2 ret;
			getMapping().toSquare(dir, ret);
			return ret;
		}
		

		Float gatherAreaPdfGMM(Vector3 wo, Float radius, 
			Vector2* componentCDFs, Vector2* componentBounds, int &topComponentCDFs, int &topComponentBounds,
			int baseCDFs, int baseBounds) const{
			// initiate sampling components
			int numNode = storedLobes;
			int pnode0 = topComponentCDFs;
			componentCDFs[topComponentCDFs].y = *(float*)&numNode;
			topComponentCDFs += numNode + 1;

			// local bound
			int numCriticalPoints = 0;
			Vector2 criticalPoints[6];
			Vector3 woLocal = localFrame.toLocal(wo);
			Float dist = woLocal.length();
			if (dist > radius){
				// project to phi plane
				Vector3 woProj = Vector3(woLocal.x, woLocal.y, 0.f);
				Float distProj = woProj.length();
				Float sqrDisTangent = distProj * distProj - radius * radius;
				// project to theta plane
				woLocal /= dist;
				Float dTheta = acos(sqrt(dist * dist - radius * radius) / dist);
				Float theta = acos(woLocal.z);
				Float theta0 = std::max(0.f, theta - dTheta);
				Float theta1 = theta + dTheta;
				if (sqrDisTangent > 0.f){
					// not covering north pole, theta-phi bounding
					Float cosdPhi = sqrt(sqrDisTangent) / distProj;
					Float dPhi = acos(cosdPhi);
					Float phi = atan2(woProj.y, woProj.x);
					Float phi0 = phi - dPhi;
					Float phi1 = phi + dPhi;			

					Float costheta0, costheta1, cosphi0, cosphi1;
					Float sintheta0, sintheta1, sinphi0, sinphi1;
					Ff::sincos(theta0, sintheta0, costheta0);
					Ff::sincos(theta1, sintheta1, costheta1);
					Ff::sincos(phi0, sinphi0, cosphi0);
					Ff::sincos(phi1, sinphi1, cosphi1);					

					criticalPoints[0] = fromPolarToSquareThetaPhi(costheta0, sintheta0, cosphi0, sinphi0);					criticalPoints[1] = fromPolarToSquareThetaPhi(costheta0, sintheta0, cosphi1, sinphi1);					criticalPoints[2] = fromPolarToSquareThetaPhi(costheta1, sintheta1, cosphi0, sinphi0);					criticalPoints[3] = fromPolarToSquareThetaPhi(costheta1, sintheta1, cosphi1, sinphi1);					numCriticalPoints = 4; 
					for (int i = 0; i < 8; i++){
						// cover the mapping changing point at 1/4 PI, add corner points
						Float deltaPhi = (-1.75f + 0.5f * (Float)i) * IMP_PI;
						if (deltaPhi > phi0 && deltaPhi < phi1){
							Float sinphi, cosphi;
							Ff::sincos(deltaPhi, sinphi, cosphi);
							criticalPoints[4] = fromPolarToSquareThetaPhi(costheta0, sintheta0, cosphi, sinphi);
							criticalPoints[5] = fromPolarToSquareThetaPhi(costheta1, sintheta1, cosphi, sinphi);
							numCriticalPoints = 6;
							break;
						}
						if (deltaPhi > phi1) break;
					}
				}
				else{ 
					// covering north pole, theta bounding
					Float costheta1, sintheta1;
					Ff::sincos(theta1, sintheta1, costheta1);
					criticalPoints[0] = fromPolarToSquareTheta(costheta1, sintheta1, 0.25f * IMP_PI);
					criticalPoints[1] = fromPolarToSquareTheta(costheta1, sintheta1, 0.75f * IMP_PI);
					criticalPoints[2] = fromPolarToSquareTheta(costheta1, sintheta1, 1.25f * IMP_PI);
					criticalPoints[3] = fromPolarToSquareTheta(costheta1, sintheta1, 1.75f * IMP_PI);
					numCriticalPoints = 4;
				}
			}
			else{
				// uniform sampling without bounding
			}

			// calculate bounding and pdf for each lobe
			Float totalProb = 0.f;
			int actualGroup = 0, actualLobe = 0;
			for (int i = 0; i < numNode; i++){
				int ptrBound = -(topComponentBounds + baseBounds);
				componentCDFs[pnode0 + i + 1].y = *(float*)&ptrBound;
				Float probLobe = lobes[actualGroup].gatherAreaPdfLobe(actualLobe, criticalPoints, numCriticalPoints, componentBounds, topComponentBounds);
				Float probi = probLobe * lobes[actualGroup].weights[actualLobe];
				componentCDFs[pnode0 + i + 1].x = probi;
				totalProb += probi;
				actualLobe++;
				if (actualLobe == TLobeType::WIDTH) {
					actualLobe = 0;
					actualGroup++;
				}
			}
			componentCDFs[pnode0].x = totalProb;
			Float invTotalProb = 1.f / totalProb;
			for (int i = 0; i < numNode; i++){
				componentCDFs[pnode0 + i + 1].x *= invTotalProb;
			}
			return totalProb;
		}

		Vector3 sampleGatherAreaGMM(Vector2 samples, Vector3 wo, Float radius, int ptrNode,
			Vector2* componentCDFs, Vector2* componentBounds) const{
			// sample lobes CDF
			Vector2 rootnode = componentCDFs[ptrNode];
			int numNode = *(int*)&rootnode.y;
			Float cdfi = 0.f;
			int chosenLobe = -1;
			Vector2 bound0, bound1;
			for (int i = 0; i < numNode; i++){
				Vector2 nodei = componentCDFs[ptrNode + 1 + i];
				Float pdfi = nodei.x;
				if (samples.x <= (cdfi + pdfi) ){
					// choose this component
					chosenLobe = i;
					int ptrBound = -*(int*)&nodei.y;					
					bound0 = componentBounds[ptrBound];
					bound1 = componentBounds[ptrBound + 1];
					samples.x = (samples.x - cdfi) / pdfi;
					break;
				}
				cdfi += pdfi;
			}

			int actualGroup = chosenLobe / TLobeType::WIDTH;
			int actualLobe = chosenLobe - actualGroup * TLobeType::WIDTH;

			const Vector2 point = lobes[actualGroup].sampleGatherArea(actualLobe, samples, bound0, bound1);
			Vector3 res;
			getMapping().fromSquare(localFrame, point, res);
			IMPORTANCE_ASSERT(res.isReal());
			return res;
		}

        IMPORTANCE_INLINE Vector3 sampleDirection(Vector2 random) const {
            const float searched = random.y;
            float cdf = 0.f, prev_cdf = 0.0f;
            int actualGroup = 0, actualLobe = 0;
            int i = 0;
            while(true) {
                prev_cdf = cdf;
                cdf += lobes[actualGroup].weights[actualLobe];
                if(cdf >= searched || i == storedLobes-1) {
                    break;
                }
                actualLobe++;
                if(actualLobe == TLobeType::WIDTH) {
                    actualLobe = 0;
                    actualGroup++;
                }
                ++i;
            }
            IMPORTANCE_ASSERT(i != storedLobes);
            //if ( i == storedLobes ) {
            //    cdf = 1.0f;
            //}
            random.y = std::min(1-FLT_EPSILON, (searched - prev_cdf)/(cdf - prev_cdf));
            IMPORTANCE_ASSERT(isReal(random.y) && random.y >= 0.f && random.y <= 1.f);

            const Vector2 point = lobes[actualGroup].sample(actualLobe, random);
            IMPORTANCE_ASSERT( lobes[actualGroup].pdf( TScalar(point.x), TScalar(point.y) )[ actualLobe ] > 0.f );
            Vector3 res;
            getMapping().fromSquare( localFrame, point, res );
            IMPORTANCE_ASSERT( res.isReal() );
            //IMPORTANCE_ASSERT( pdf( res ) > 0.f );
            return res;
        }

        IMPORTANCE_INLINE float pdf(const Vector3 dir) const {
            Vector2 x;            
            if ( getMapping().toSquare( localFrame, dir, x ) == false ) {
                return 0.0f;
            }            
            float res = 0.f;
            const int cnt = (storedLobes+TLobeType::WIDTH-1) / TLobeType::WIDTH;
            for (int i = 0; i < cnt; ++i ) {
                res += TScalar::dot( lobes[ i ].weights, lobes[i].pdf(TScalar(x.x), TScalar(x.y)) );
            }
            return res / TMapping::jacobian();
        }

        static IMPORTANCE_INLINE TMapping getMapping() {
            return TMapping();
        }
    };

    typedef GaussianSamplingMix<MAX_FITTED_LOBES, GaussianSamplingLobes<Float4>, SquareHemisphereMapping> LerpGaussian;


	template<typename TScalar>
    class GaussianMultiLobe : public GaussianSamplingLobes<TScalar> {
    public:
		
        // Online stepwise E-M        
        TScalar statsW, statsGamma_w;
		TScalar x1, x2, x1sqr, x2sqr, x1x2;
    
        IMPORTANCE_INLINE GaussianMultiLobe() {
			cov.identity( TScalar( 0.0f ) );
			normalization = mean.x = mean.y = weights = 
			statsW = statsGamma_w = TScalar(0.0f);
			x1 = x2 = x1sqr = x2sqr = x1x2 = TScalar( 0.f );
        }

        std::string toString(int index) const {
            std::stringstream str;
            str << "mean = " << Vector2(mean.x[index],mean.y[index]).toString() << " covariance matrix: [ " << cov.m[ 0 ][index] << ", " << cov.m[ 1 ][index] << " | " << cov.m[ 1 ][index] << ", " << cov.m[ 2 ][index] << " ]";
            return str.str();
        }
                
		// Changes cov_** -> should be called only once !
		IMPORTANCE_INLINE void init() {
            const TScalar det_cov = cov.determinant();
			IMPORTANCE_ASSERT( det_cov.isReal() );
			cov.invert();
			normalization = TScalar(1.0f) / ( TScalar( 2.0f * IMP_PI ) * det_cov.abs().sqrt() ); 
			IMPORTANCE_ASSERT( normalization.isReal() );
			IMPORTANCE_ASSERT( cov.m[ 0 ].isReal() );
			IMPORTANCE_ASSERT( cov.m[ 1 ].isReal() );
			IMPORTANCE_ASSERT( cov.m[ 2 ].isReal() );
        }    
        
		IMPORTANCE_INLINE float cachingStdDev(const int index) const {
            return 1.0f / sqrt(sqrt(cov.determinant()[ index ]));
		}

        IMPORTANCE_INLINE float sigmaInvertedEigenvalue(const int index) const {
            const float trace = cov.trace()[index];
            const float sqrtTerm = trace*trace - 4*cov.determinant()[index];
            //IMPORTANCE_ASSERT(sqrtTerm > -1e-3f);
            const float result = (trace + sqrt(std::max(sqrtTerm, 0.f))) /2;
            IMPORTANCE_ASSERT(isReal(result));
            return result;
        }

    };

	template<int TLobesCount, int TVectorWidth, typename TScalar, typename TMapping>
    class GaussianMultiMix : public GaussianSamplingMix<TLobesCount, GaussianMultiLobe<TScalar>, TMapping> {
		class WeightComparator {
            WeightComparator& operator=(const WeightComparator&);
        public:
            WeightComparator( const GaussianMultiMix & gmm ) : m_gmm( gmm ) {}
            bool operator() ( int i,int j ) { return ( m_gmm.getWeight( i ) > m_gmm.getWeight( j ) ); }
        private:
            const GaussianMultiMix & m_gmm;
        };
	public:
		static const TScalar EPSILON;
		static const TScalar MIN_DETERMINANT;        
        static const Float	 MINIMAL_ERROR;		
        static const size_t MINIMAL_DATASIZE = 5;		
        static const TScalar PRIOR_A;
        static const TScalar PRIOR_B;
        static const TScalar PRIOR_V;
		typedef GaussianMultiLobe<TScalar> MultiLobe;

		/// Online EM statistics        
        unsigned int k;             // number of points
        Float statsPointWeight;     // total weight of all data points        
		bool converged;				// has the distribution converged ?

#ifdef LIBIMP_STATS
		EmInfo info;
#define EMINFO(gmm) &(gmm).info
#else
#define EMINFO(x) NULL
#endif

        IMPORTANCE_INLINE GaussianMultiMix( int lobes = 0 ) {
            storedLobes			= lobes;
			k					= 0;
            statsPointWeight	= 0.f;
			converged			= false;            
        }


IMPORTANCE_INLINE const EmInfo * getInfo() const { 
#ifdef LIBIMP_STATS
    return &info; 
#else
    return NULL;
#endif
}

#ifdef LIBIMP_STATS
virtual void release() { info.release(); }
#endif

        std::string toString() const {
			std::vector<int> indices;
            for ( int i = 0; i < nComponents(); ++i ) { indices.push_back( i ); }
            std::sort( indices.begin(), indices.end(), WeightComparator( *this ) );
			
            std::stringstream str;
            str << "GaussianMultiMix [ " << std::endl;

			for ( int h = 0; h < storedLobes; ++h ) {
				const div_t tmp = div(h, TVectorWidth);
                int index = indices[ h ];
                str << "\tC[ " << index << " ]: " << "w: "<< getWeight( index ) << " " << lobes[tmp.quot].toString(tmp.rem) << std::endl;
            }

            str << std::endl 
                << "k = " << k << std::endl
                << "statsPointWeight = " << statsPointWeight << std::endl
                << "]" << std::endl;

            return str.str();
		}

		/* Check weights - have to sum to one, only for debugging */
        IMPORTANCE_INLINE bool check( Float & error, Float & weightsSum ) {        
            weightsSum = 0.0f;   
			const int cnt = (storedLobes+TVectorWidth-1) / TVectorWidth;
			static const TScalar v1( 1.0f );
            for (int i = 0; i < cnt; ++i ) {
                weightsSum += TScalar::dot( lobes[ i ].weights, v1 );
            }            
            error = std::abs( 1.0f - weightsSum );
            return error < 1e-3;
        }

        // podle http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence#KL_divergence_for_Normal_Distributions
        static IMPORTANCE_INLINE Float klDist(const MultiLobe& lobe0, const int index0, const MultiLobe& lobe1, const int index1, Vector2 meanDiff ) {
			const SPDMatrix2x2< Float> inv0(lobe0.cov.m[0][ index0 ], lobe0.cov.m[1][ index0 ], lobe0.cov.m[2][ index0 ] );
			const SPDMatrix2x2< Float> inv1(lobe1.cov.m[0][ index1 ], lobe1.cov.m[1][ index1 ], lobe1.cov.m[2][ index1 ] );
            const SPDMatrix2x2< Float > cov0 = inv0.getInverse();
            const SPDMatrix2x2< Float > cov1 = inv1.getInverse();
            
            // prvky ktery potrebuju z inv1*cov0
            const Float mult0 = inv1.m[0]*cov0.m[0] + inv1.m[1]*cov0.m[1];
            const Float mult1 = inv1.m[1]*cov0.m[1] + inv1.m[2]*cov0.m[2];
            const Float trace = mult0+mult1;

            
            const Float norm = inv1.quadraticForm(meanDiff);
            const Float result = (trace + norm - 2 - Ff::log(cov0.determinant()/cov1.determinant()))/2.f;
            IMPORTANCE_ASSERT(isReal(result));
            return result;
        }
        static IMPORTANCE_INLINE Float klDist(const MultiLobe& lobe0, const int index0, const MultiLobe& lobe1, const int index1 ) {
            const Vector2 meanDiff( lobe1.mean.x[ index1 ] - lobe0.mean.x[ index0 ], lobe1.mean.y[ index1 ] - lobe0.mean.y[ index0 ] );
            return klDist(lobe0, index0, lobe1, index1, meanDiff);
        }
        
        IMPORTANCE_INLINE static Float lobeDistance(const GaussianMultiMix& a, const int indexA, const GaussianMultiMix& b, const int indexB, const bool forMixture) {
			const div_t tmpA = div(indexA, TVectorWidth);
			const div_t tmpB = div(indexB, TVectorWidth);
            const Float wDiff = a.lobes[tmpA.quot].weights[tmpA.rem] - b.lobes[tmpB.quot].weights[tmpB.rem];

            const Vector2 mean0 = Vector2(a.lobes[tmpA.quot].mean.x[tmpA.rem], a.lobes[tmpA.quot].mean.y[tmpA.rem]);
            const Vector2 mean1 = Vector2(b.lobes[tmpB.quot].mean.x[tmpB.rem], b.lobes[tmpB.quot].mean.y[tmpB.rem]);

            Frame rotMatrix = a.localFrame * b.localFrame.getInverse();
            Vector2 newMean1;
            Vector3 dir;
            getMapping().fromSquare( mean1, dir );                
            getMapping().toSquare(rotMatrix.toLocal( dir ), newMean1 );

            const Vector2 meanDiff = mean0-newMean1;

            if(forMixture) {
                const Float kl1 = klDist(a.lobes[tmpA.quot], tmpA.rem, b.lobes[tmpB.quot], tmpB.rem, meanDiff);
                IMPORTANCE_ASSERT(isReal(kl1));
                return kl1;
            } else {
                const Float kl1 = klDist(a.lobes[tmpA.quot], tmpA.rem, b.lobes[tmpB.quot], tmpB.rem, meanDiff);
                const Float kl2 = klDist(b.lobes[tmpB.quot], tmpB.rem, a.lobes[tmpA.quot], tmpA.rem, meanDiff);
                const Float result = kl1 + kl2 + 2*abs(wDiff);
                IMPORTANCE_ASSERT(isReal(result));
                return result;
            }
        }

        // http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.148.2502&rep=rep1&type=pdf eq 20
        IMPORTANCE_INLINE static Float klDistance(const GaussianMultiMix& f, const GaussianMultiMix& g) {

            Float result = 0;
            for(int a = 0; a < f.nComponents(); ++a) {
				const div_t tmpA = div(a, TVectorWidth);
                Float top = 0;
                for(int it_a = 0; it_a < f.nComponents(); ++it_a) {
                    const div_t _a = div(it_a, TVectorWidth);
                    top += f.lobes[_a.quot].weights[_a.rem] * Ff::exp( -lobeDistance(f, a, f, it_a, true));
                }
                Float bottom = 0;
				for(int it_b = 0; it_b < g.nComponents(); ++it_b) {
                    const div_t b = div(it_b, TVectorWidth);
                    bottom += g.lobes[ b.quot ].weights[b.rem] * Ff::exp(-lobeDistance(f, a, g, it_b, true));
                }
                result += f.lobes[tmpA.quot].weights[tmpA.rem] * log(top/std::max(1e-6f, bottom));
				IMPORTANCE_ASSERT(isReal(result));
            }
            IMPORTANCE_ASSERT(isReal(result));
            return result;
        }

		static void mixtureLerp(const GaussianMultiMix &a, const GaussianMultiMix &b, const unsigned char * pairing, const float amountB, LerpGaussian& out) {
            IMPORTANCE_ASSERT(isReal(amountB));
            
            out.localFrame = a.localFrame;
            out.storedLobes = a.nComponents();
            const int activeSseUnits = (out.storedLobes+LerpGaussian::WIDTH-1)/LerpGaussian::WIDTH;

            ImStaticArray<Float4, MAX_FITTED_LOBES/LerpGaussian::WIDTH> meanX, meanY;

            if(dot(a.localFrame.n, b.localFrame.n) > 0.997f && dot(a.localFrame.s, b.localFrame.s) > 0.997f) {
                for(int i = 0; i < activeSseUnits; ++i) {
                    meanX[i] = b.lobes[i].mean.x;
                    meanY[i] = b.lobes[i].mean.y;
                }
            } else {
                Frame rotMatrix = a.localFrame * b.localFrame.getInverse();
                for(int i = 0; i < activeSseUnits; ++i) {
                    Vector3 dir;
                    Vector2 point;
                    getMapping().fromSquare( Vector2( b.lobes[i].mean.x, b.lobes[i].mean.y ), dir );
                    getMapping().toSquare( rotMatrix.toLocal(dir), point );
                    meanX[i] = point.x; 
                    meanY[i] = point.y;
                    //IMPORTANCE_ASSERT(outLobes[i].mean.isReal()));
                }
            }
            

			// Handle pairings:
            int ptr = 0;
			for ( int i = 0; i < activeSseUnits; ++i ) {
                for(int j = 0; j < 4; ++j) {

				    const div_t tmpB = div(pairing[ptr++], LerpGaussian::WIDTH);

				    out.lobes[ i ].mean.x    [ j ] = meanX  [ tmpB.quot ]           [ tmpB.rem ];
				    out.lobes[ i ].mean.y    [ j ] = meanY  [ tmpB.quot ]           [ tmpB.rem ];
				    out.lobes[ i ].cov.m[ 0 ][ j ] = b.lobes[ tmpB.quot ].cov.m[ 0 ][ tmpB.rem ];
				    out.lobes[ i ].cov.m[ 1 ][ j ] = b.lobes[ tmpB.quot ].cov.m[ 1 ][ tmpB.rem ];
				    out.lobes[ i ].cov.m[ 2 ][ j ] = b.lobes[ tmpB.quot ].cov.m[ 2 ][ tmpB.rem ];
				    out.lobes[ i ].weights   [ j ] = b.lobes[ tmpB.quot ].weights   [ tmpB.rem ];
                }
			}

			// Lerp
			for ( int i = 0; i < activeSseUnits; ++i ) {
				out.lobes[ i ].mean.x = lerp( a.lobes[ i ].mean.x, out.lobes[ i ].mean.x, amountB );
				out.lobes[ i ].mean.y = lerp( a.lobes[ i ].mean.y, out.lobes[ i ].mean.y, amountB );
                out.lobes[ i ].weights = lerp( a.lobes[ i ].weights, out.lobes[ i ].weights, amountB );
                //out.lobes[ i ].cov = SPDMatrix2x2< TScalar >::lerp( a.lobes[i].cov.getInverse(), out.lobes[ i ].cov.getInverse(), amountB ).getInverse();
				out.lobes[ i ].cov = SPDMatrix2x2< TScalar >::lerp( a.lobes[i].cov, out.lobes[ i ].cov, amountB );
				out.lobes[ i ].lightweightInit();
			}
#ifdef LIBIMP_DEBUG
            float sumW = 0;
            for(int i = 0; i < activeSseUnits; ++i) {
                for(int j = 0; j < 4; ++j) {
                    sumW += out.lobes[i].weights[j];
                }
            }
            IMPORTANCE_ASSERT(abs(sumW-1) < 1e-3f);
#endif

		}
        IMPORTANCE_INLINE int nComponents() const {
            return storedLobes;
        }

        IMPORTANCE_INLINE float getWeight(const int index) const {
            const div_t tmp = div(index, TVectorWidth);
            return lobes[tmp.quot].weights[tmp.rem];
        }

        
        

        IMPORTANCE_INLINE float sigmaInvertedEigenvalue(const int i) const
        {
            const div_t tmp = div(i, TVectorWidth);
            return lobes[tmp.quot].sigmaInvertedEigenvalue(tmp.rem);
        }

		IMPORTANCE_INLINE const Vector3 lobeLocalDir(const int i) const
		{
            const div_t tmp = div(i, TVectorWidth);
            Vector3 dir;
            getMapping().fromSquare( Vector2( lobes[tmp.quot].mean.x[tmp.rem], lobes[tmp.quot].mean.y[tmp.rem] ), dir );
			return dir;
		}

		IMPORTANCE_INLINE void randomInit( Importance::Random & rnd ) {
			k					= 0;
            statsPointWeight	= 0.f; 
			for ( int i = 0; i < storedLobes; ++i ) {
				const div_t tmp = div(i, TVectorWidth);
				lobes[ tmp.quot ].mean.x[ tmp.rem ] = Importance::nextFloat(rnd);
				lobes[ tmp.quot ].mean.y[ tmp.rem ] = Importance::nextFloat(rnd);
			}
			const int cnt = storedLobes / TVectorWidth;
			const TScalar weight = TScalar(1.0f / storedLobes);
			for ( int i = 0; i < cnt; ++i ) {
				lobes[ i ].cov.identity( TScalar(0.1f) );
				lobes[ i ].weights = weight;
				lobes[ i ].init();
			}
		}

        IMPORTANCE_INLINE void randomInitMcLachlan( Importance::Random & rnd, const ArrayGaussianSSE< TScalar, TVectorWidth, GaussianMultiMix> & samples ) {
            IMPORTANCE_ASSERT( samples.size() > 2 );
            QuadVector2 mean( TScalar( 0.f ) );            
            for ( size_t i = 0; i < samples.size(); ++i ) {
                mean = mean + samples[ i ].getPoint();
            }
            mean = mean / TScalar( (float) samples.size() );
            
            SPDMatrix2x2<TScalar> cov;                        
            for ( size_t i = 0; i < samples.size(); ++i ) {        
                QuadVector2 p = samples[ i ].getPoint() - mean;        
                cov += SPDMatrix2x2<TScalar>( p );
            }
            cov = cov / TScalar( (float) samples.size() );

            cov.invert();
            GaussianSamplingLobes<TScalar> normDist;
            normDist.cov  = cov;
            normDist.mean = mean;       
            normDist.lightweightInit();

            for ( int i = 0; i < storedLobes; ++i ) {
                Vector2 sm( -1.f, -1.f );
                int j = 0;
                while ( j < 4 && ( sm.x < 0.0f || sm.x > 1.0f || sm.y < 0.0f || sm.y > 1.0f ) ) {
                    sm = normDist.sampleSquare( 0, Vector2( nextFloat( rnd ), nextFloat( rnd ) ) );
                    j ++;
                }
                if ( j > 3 ) {
                    IMPORTANCE_ASSERT( true );
                    sm.x = clamp( sm.x, 0.f, 1.f );
                    sm.y = clamp( sm.y, 0.f, 1.f );
                }
                const div_t tmp = div( i, TVectorWidth );
                lobes[ tmp.quot ].mean.x[ tmp.rem ] = sm.x;
                lobes[ tmp.quot ].mean.y[ tmp.rem ] = sm.y;                                
            }
        
            const TScalar weight = TScalar(1.0f / storedLobes);
            for ( int i = 0; i < storedLobes / TVectorWidth; ++i ) {
                lobes[ i ].cov      = cov;
                lobes[ i ].weights  = weight;
                lobes[ i ].init();
            }
        }

		IMPORTANCE_INLINE void initAsUninformed() {
			storedLobes = 4;
			lobes[ 0 ].mean.x = TScalar(0.5f);
			lobes[ 0 ].mean.y = TScalar(0.5f);
			lobes[ 0 ].cov.identity( TScalar(1.0f) );
			lobes[ 0 ].weights = TScalar(0.25f);
			lobes[ 0 ].avgDistance = TScalar(0.0f);
			lobes[ 0 ].init();            
		}

        //////////////////////////////////////////////////////////////////////////
        // Caching stuff
        //////////////////////////////////////////////////////////////////////////
        static IMPORTANCE_INLINE Float validityRadius(const GaussianMultiMix& dist, const Config& config, const float pixel2world = 0.f ) {
            Float result = 0;
            Float weights = 0.f;
            Float avg = 0;
            Float minimum = 0;
            Float totalDistance = dist.m_cacheStats.avgParticleDistance;

            for(int i = 0; i < dist.nComponents(); i++) {
                const float eigenVal = dist.sigmaInvertedEigenvalue(i);
                //TODO: rename minProjectedIntensity to KLDivergence threshold
                const float maxProjectedDiff = sqrt(config.cache.klDivThres / eigenVal);

                const Vector2 tmpPoint(std::min(1.f, maxProjectedDiff+0.5f), 0.5f);
                Vector3 dir;                
                getMapping().fromSquare( tmpPoint, dir );               
                const float cosMaxAngle = clamp(dir.z, 1e-2f, 0.99999f);
                const float maxAngle = acos(cosMaxAngle);

                const float w = dist.getWeight(i);
                const float maxDifference = tan(maxAngle);
                
                result += w/(maxDifference);
                weights += dist.getWeight(i);
                avg += w*maxDifference;
                minimum = std::min(maxDifference, minimum);
            }

            avg *= totalDistance;
            minimum *= totalDistance;
            result = totalDistance/result;

            IMPORTANCE_ASSERT(abs(weights-1.f) < 1e-4f);
            IMPORTANCE_ASSERT(result <= avg*1.001f);
            IMPORTANCE_ASSERT(result >= minimum);
            //return 1/result;
            //return avg;
            return result;
        }
    };

}

#pragma warning(pop)