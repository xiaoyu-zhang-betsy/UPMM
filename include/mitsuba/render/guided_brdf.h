/*
    This file is part of a demo implementation of an importance sampling technique
    described in the "On-line Learning of Parametric Mixture Models for Light Transport Simulation"
    (SIGGRAPH 2014) paper.
    The implementation is based on Mitsuba, a physically based rendering system.

    Copyright (c) 2014 by Jiri Vorba, Ondrej Karlik, Martin Sik.
    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/libImpUtils.h>

MTS_NAMESPACE_BEGIN

/** Easy sampling from MIS combination of bsdf and trained radiace/importance distribution */
class GuidedBRDF {
public:
    /** Constructs guided sampler */
	GuidedBRDF(
		mitsuba::Intersection       its,
		Importance::Sampler*        importanceSampler,
		Float                       bsdfSamplingProbability,
		bool valid = true)
		: m_its(its), m_importanceSampler(importanceSampler),
		m_bsdfSamplingProbability(bsdfSamplingProbability), m_impDistrib(NULL), m_eta(-1.f) {

		if (!valid) return;

		/* Prepare distribution if guided path-tracing is set on*/
		m_bsdf = its.getBSDF();
		if (m_importanceSampler != NULL) {
			Importance::Hit hit = itsToHit(its);
			m_impDistrib = m_importanceSampler->getDistribution(hit, m_distBuffer);
#ifdef LIBIMP_STATS
			if (m_impDistrib != NULL) {
				/// somehow prevent storing photons from EM if you don't need it and statistics are on
				Importance::EmInfo * info = const_cast<Importance::EmInfo *>(m_impDistrib->getInfo());
				if (info != NULL) {
					info->particles.clear();
				}
			}
#endif
		}
		m_pureSpecularBSDF = (m_bsdf->getType() & BSDF::EAll & ~BSDF::EDelta) == 0;
	}

    GuidedBRDF( 
        mitsuba::Intersection       its, 
        const RayDifferential&      ray, 
        Importance::Sampler*        importanceSampler,                     
        Float                       bsdfSamplingProbability ) 
        : m_its( its ), m_importanceSampler( importanceSampler),
        m_bsdfSamplingProbability( bsdfSamplingProbability ), m_impDistrib( NULL ), m_eta( -1.f ) {

            /* Prepare distribution if guided path-tracing is set on*/
            m_bsdf = its.getBSDF(ray);
            if ( m_importanceSampler != NULL ) {
                Importance::Hit hit = itsToHit( its );
                m_impDistrib = m_importanceSampler->getDistribution( hit, m_distBuffer );
#ifdef LIBIMP_STATS
                if ( m_impDistrib != NULL ) {
                    /// somehow prevent storing photons from EM if you don't need it and statistics are on
                    Importance::EmInfo * info = const_cast<Importance::EmInfo *>( m_impDistrib->getInfo() );                
                    if ( info != NULL ) {                    
                        info->particles.clear();
                    }
                }
#endif
            }
            m_pureSpecularBSDF = ( m_bsdf->getType() & BSDF::EAll & ~BSDF::EDelta ) == 0;
    }

    /// Returns true, if bsdf is purely specular (i.e. only delta transmission and/or reflection)
    bool isBSDFPureSpecular() const {
        return m_pureSpecularBSDF;
    }

    /// Returns false, if we failed to construct distribution
    bool isValid() const {
        return m_importanceSampler == NULL || m_impDistrib != NULL;
    }

    /// Evaluates bsdf value at given an outgoing direction (direction is in world space!)
    Spectrum eval( const Vector & direction ) const {
        /* Evaluate BSDF * cos(theta) */
        BSDFSamplingRecord bRec(m_its, m_its.toLocal(direction));
        return m_bsdf->eval(bRec);
    }

    /// Returns a pdf of combined bsdf & guided sampling given an outgoing direction (direction is in world space!)
    float pdf( const Vector & direction) const {
        /* Calculate prob. of having generated that direction using BSDF sampling */
        BSDFSamplingRecord bRec(m_its, m_its.toLocal(direction));
        Float bsdfPdf = m_bsdf->pdf(bRec);

        Float iiPdf = 0.0f;
        /// count in IIS strategy
        if ( m_impDistrib != NULL && !m_pureSpecularBSDF ) {
            m_impDistrib->pdfs( &IMP_VECTOR3( direction ), &iiPdf, 1 );                    
            iiPdf *= (1 - m_bsdfSamplingProbability);
            bsdfPdf *= m_bsdfSamplingProbability;
        }
        return iiPdf + bsdfPdf;
    }

    /// Samples a new outgoing direction, returns bsdf value divided by pdf, the new direction, pdf, sampled type (if bsdf was used, otherwise returns 0), bsdf eta
    /// Uses MIS with balance heuristic to combine sampling from bsdf and trained distribution
    Spectrum sample( Vector & direction, Float & pdf, Sampler *sampler) {
        Spectrum albedo;
        return sample( direction, pdf, sampler, albedo);
    }

    /// Samples a new outgoing direction, returns bsdf value divided by pdf, the new direction, pdf, sampled type (if bsdf was used, otherwise returns 0), bsdf eta and albedo
    /// Uses MIS with balance heuristic to combine sampling from bsdf and trained distribution
    Spectrum sample( Vector & direction, Float & pdf, Sampler *sampler, Spectrum & albedo ) {
        /* Determine whether to sample according to BSDS or according to illumination */
        bool sampleOnlyBSDF = m_pureSpecularBSDF || m_impDistrib == NULL;
        bool useBSDF = sampleOnlyBSDF || sampler->next1D() < m_bsdfSamplingProbability;

        if ( useBSDF ) {
            /* Sample BSDF * cos(theta) */
            BSDFSamplingRecord bRec(m_its, sampler, ERadiance);
            Spectrum out = albedo = m_bsdf->sample(bRec, pdf, sampler->next2D());
            if ( out.isZero() ) {
                return Spectrum( 0.f );
            }
            out *= pdf; /* Already pre-multiplied by inverse of mainPdf */
            m_lastSampledComponent = bRec.sampledType;
            m_eta = bRec.eta;
            direction = m_its.toWorld( bRec.wo );
            // Sampling from distribution is possible?
            if ( !sampleOnlyBSDF ) {
                pdf *= m_bsdfSamplingProbability; // Multiply by probability of selecting bsdf as sampling strategy
                if ((!(bRec.sampledType & BSDF::EDelta)) ) {// Not dirac surface				
                    Float iiPdf;
                    m_impDistrib->pdfs( &IMP_VECTOR3( direction ), &iiPdf, 1 );
                    pdf += iiPdf * ( 1 - m_bsdfSamplingProbability ); // Combined pdf
                }
            }
            if ( pdf != 0.0f ) {
                return out / pdf;
            } else {
                return Spectrum(0.0f);
            }
        } else {
            // Sample the distribution
            Importance::Vector3 res;
            Point2 s1 = sampler->next2D();                
            Importance::Vector2 samples( s1.x, s1.y);
            m_impDistrib->sampleDirections( &samples, &res, &pdf, 1 );
            pdf *= ( 1.0f - m_bsdfSamplingProbability );
            direction = MTS_VECTOR( res );
            /// Evaluate BRDF in the new direction * cos(theta) 
            BSDFSamplingRecord bRec( m_its, m_its.toLocal( direction ) );
            m_eta = m_bsdf->getEta();
            m_lastSampledComponent = 0;
            Float bsdfPdf = m_bsdf->pdf( bRec );
            pdf += bsdfPdf * m_bsdfSamplingProbability;
            if ( pdf != 0.0f ) {
                Spectrum out = albedo = m_bsdf->eval( bRec );
                if ( bsdfPdf != 0.0f ) {
                    albedo /= bsdfPdf;
                } else {
                    albedo = Spectrum(0.0f);
                }
                return  out / pdf;
            } else {
                albedo = Spectrum(0.0f);
                return albedo;
            }
        }
    }

	Spectrum sampleGMM(Vector & direction, Float & pdf, Sampler *sampler) {
		// Sample the distribution
		Importance::Vector3 res;
		Point2 s1 = sampler->next2D();
		Importance::Vector2 samples(s1.x, s1.y);
		m_impDistrib->sampleDirections(&samples, &res, &pdf, 1);
		pdf *= (1.0f - m_bsdfSamplingProbability);
		direction = MTS_VECTOR(res);
		/// Evaluate BRDF in the new direction * cos(theta) 
		BSDFSamplingRecord bRec(m_its, m_its.toLocal(direction));
		m_eta = m_bsdf->getEta();
		m_lastSampledComponent = 0;
		Float bsdfPdf = m_bsdf->pdf(bRec);
		pdf += bsdfPdf * m_bsdfSamplingProbability;
		if (pdf != 0.0f) {
			Spectrum out = m_bsdf->eval(bRec);
			return  out / pdf;
		}
		else {
			return Spectrum(0.0f);
		}
	}

    /// Type of component sampled using last bsdf sampling
    unsigned int getLastSampledComponent() const {
        return m_lastSampledComponent;
    }

    /// Eta at given intersection
    float getEta() const {
        SAssert( m_eta > -1.f );
        return m_eta;
    }

    /// Destructor
    ~GuidedBRDF() {
        if ( m_impDistrib != NULL ) {
            m_impDistrib->release();
            m_impDistrib = NULL;
        }
    }
public:
    /// Intersection from which we want to generate next ray
    mitsuba::Intersection & m_its;
    /// BSDF at intersection
    const BSDF *m_bsdf;
    /// Radiance/importance distribution sampler
    Importance::Sampler * m_importanceSampler;
    /// Buffer where the current distribution is allocated to avoid expensive memory re-allocation
    Importance::ResultBuffer m_distBuffer;
    /// The distribution will be allocated into the distBuffer so do not dispose this memory by delete!!
    Importance::Distribution * m_impDistrib;
    /// Probability to sample the next ray from bsdf
    Float m_bsdfSamplingProbability;
    /// Non-dirac bsdf?
    bool m_pureSpecularBSDF;
    /// Type of component sampled using last bsdf sampling
    unsigned int m_lastSampledComponent;
    /// Eta at given intersection
    float m_eta;        
};

MTS_NAMESPACE_END