/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2012 by Wenzel Jakob and others.

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

#include <mitsuba/render/guided_gatherproc.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/guided_brdf.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief This work result implementation stores a sequence of photons, which can be
 * sent over the wire as needed.
 *
 * It is used to implement parallel networked photon tracing passes.
 */
class GuidedPhotonVector : public WorkResult {
public:
	GuidedPhotonVector() { }

	inline void nextParticle() {
		m_particleIndices.push_back((uint32_t) m_photons.size());
	}

	inline void put(const Photon &p) {
		m_photons.push_back(p);
	}

	inline size_t size() const {
		return m_photons.size();
	}

	inline size_t getParticleCount() const {
		return m_particleIndices.size()-1;
	}

	inline size_t getParticleIndex(size_t idx) const {
		return m_particleIndices.at(idx);
	}

	inline void clear() {
		m_photons.clear();
		m_particleIndices.clear();
		m_bgParticles.m_data.clear();
	}

	inline const Photon &operator[](size_t index) const {
		return m_photons[index];
	}

	inline std::vector<Photon>::const_iterator bgParticlesBegin() const {
		return m_bgParticles.m_data.begin();
	}

	inline const std::vector<Photon>::const_iterator bgParticlesEnd() const {
		return m_bgParticles.m_data.end();
	}

	inline void putBgParticle( const Photon & particle ) {
		m_bgParticles.m_data.push_back( particle );
	}

	void load(Stream *stream) {
		throw std::runtime_error("Serialization is not supported!");
	}

	void save(Stream *stream) const {
		throw std::runtime_error("Serialization is not supported!");
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "GuidedPhotonVector[size=" << m_photons.size() << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
public:
    ParticleStats m_particleStats;
protected:
	// Virtual destructor
	virtual ~GuidedPhotonVector() { }
private:
	std::vector<Photon> m_photons;
	std::vector<uint32_t> m_particleIndices;
	BgParticlesVec m_bgParticles;    
};

/**
 * This class does the actual photon tracing work
 */
class GuidedGatherPhotonWorker : public ParticleTracer {
public:
	GuidedGatherPhotonWorker(
        GuidedGatherPhotonProcess::EGatherType type, 
        size_t granularity,
		int maxDepth, 
        int rrDepth, 
        const WeightWindow & ww,
        Importance::Sampler * importance,		
        const GuidingConfig & cfg,
        GuidedGatherPhotonProcess::ETraceDirection traceDirection) 
        : ParticleTracer(maxDepth, rrDepth, false),
		m_type(type), 
		m_traceDirection( traceDirection ), 
		m_importance( importance ), 
		m_granularity(granularity), 		
        m_ww( ww ),
        m_cfg( cfg ) {}

	GuidedGatherPhotonWorker(Stream *stream, InstanceManager *manager)
	 : ParticleTracer(stream, manager) {
		throw std::runtime_error("Serialization is not supported!");
	}

	ref<WorkProcessor> clone() const {
		return new GuidedGatherPhotonWorker(m_type, m_granularity, m_maxDepth,
			m_rrDepth, m_ww, m_importance, m_cfg, m_traceDirection);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		throw std::runtime_error("Serialization is not supported!");
	}

	ref<WorkResult> createWorkResult() const {
		return new GuidedPhotonVector();
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult,
		const bool &stop) {
		m_workResult = static_cast<GuidedPhotonVector *>(workResult);
		m_workResult->clear();

		const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);
		MediumSamplingRecord mRec;
		Intersection its;


		ref<Sensor> sensor    = m_scene->getSensor();
		bool needsTimeSample  = sensor->needsTimeSample();
		PositionSamplingRecord pRec(sensor->getShutterOpen()
			+ 0.5f * sensor->getShutterOpenTime());
		Ray ray;

		m_sampler->generate(Point2i(0));

		for (size_t index = range->getRangeStart(); index <= range->getRangeEnd() && !stop; ++index) {
			m_sampler->setSampleIndex(index);

			Float timeSample = 0.5f;
			/* Sample an emission */
			if ( needsTimeSample ) {
				timeSample = m_sampler->next1D();
				pRec.time = sensor->sampleTime(timeSample);
			}

			const Emitter *emitter = NULL;
			const Medium *medium;

			Spectrum power;
			Ray ray;
			Point2 pixel;

			if ( m_traceDirection == GuidedGatherPhotonProcess::EFromEmitters ) {                
				/* Sample both components together, which is potentially
					faster / uses a better sampling strategy */
				power = m_scene->sampleEmitterRay(ray, emitter,
					m_sampler->next2D(), m_sampler->next2D(), pRec.time);
				medium = emitter->getMedium();
                if ( power.isZero() ) {
                    continue;
                }

				handleNewParticle();
            } else { // Tracing from sensor
				pixel = m_sampler->next2D();
				pixel.x *= sensor->getFilm()->getCropSize().x;
				pixel.y *= sensor->getFilm()->getCropSize().y;
				pixel += Point2(sensor->getFilm()->getCropOffset());
				Point2 apertureSample = Point2(0.5f);
				if (sensor->needsApertureSample())
					apertureSample = m_sampler->next2D();
				power = sensor->sampleRay(ray, pixel, apertureSample, timeSample );
				medium = sensor->getMedium();
				handleNewParticle();
			}

			int depth               = 1, 
                nullInteractions    = 0;
			bool delta              = false;
            Float eta               = 1.f;


            while( zasobnik not empty ) {}

			Spectrum throughput(1.0f); // unitless path throughput (used for russian roulette)
			while (!throughput.isZero() && (depth <= m_maxDepth || m_maxDepth < 0)) {
				m_scene->rayIntersect(ray, its);
				Spectrum albedo;
				Spectrum bsdfWeight;


                if ( split ) {
                    ray, its, depth, medium, throughput, power, nullInteractions, delta, eta
                    uloz na zasobnik aktualni cestu podle poctu splitu
                }


				/* ==================================================================== */
				/*                 Radiative Transfer Equation sampling                 */
				/* ==================================================================== */
				if (medium && medium->sampleDistance(Ray(ray, 0, its.t), mRec, m_sampler)) {
					/* Sample the integral
					  \int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
					*/

					throughput *= mRec.sigmaS * mRec.transmittance / mRec.pdfSuccess;

					/* Forward the medium scattering event to the attached handler */
					handleMediumInteraction(depth, nullInteractions,
							delta, mRec, medium, -ray.d, throughput*power);

					PhaseFunctionSamplingRecord pRec(mRec, -ray.d, EImportance);

					throughput *= medium->getPhaseFunction()->sample(pRec, m_sampler);
					delta = false;

					ray = Ray(mRec.p, pRec.wo, ray.time);
					ray.mint = 0;
				} else if (its.t == std::numeric_limits<Float>::infinity()) {
					/* There is no surface in this direction */
					if ( m_scene->hasEnvironmentEmitter() && m_traceDirection == GuidedGatherPhotonProcess::EFromSensor ) {
						if (medium)
							throughput *= mRec.transmittance / mRec.pdfFailure;
						handleBackgroundInteraction( depth, ray, throughput * power );
					}
					break;
				} else {
					/* Sample
						tau(x, y) (Surface integral). This happens with probability mRec.pdfFailure
						Account for this and multiply by the proper per-color-channel transmittance.
					*/
					if ( medium ) {
						throughput *= mRec.transmittance / mRec.pdfFailure;
                    }
					
					Float distance = (its.p - ray.o).length();
					
					/* Forward the surface scattering event to the attached handler */
					handleSurfaceInteraction(depth, nullInteractions, delta, its, medium, throughput*power, distance);

					/* Prepare importance & bsdf sampling distribution */
					RayDifferential rayd( ray );
					GuidedBRDF gsampler( its, rayd, m_importance, m_cfg.m_mitsuba.bsdfSamplingProbability );

					/* Sample ray direction */
					Vector wo;
					Float pdf;
					Spectrum albedo;
					bsdfWeight = gsampler.sample(wo, pdf, m_sampler, albedo);

					if ( bsdfWeight.isZero() || pdf == 0.0f ) {
						break;
                    }

					/* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
					Vector wi = -ray.d;
					//Float wiDotGeoN = dot(its.geoFrame.n, wi);
					Float woDotGeoN = dot(its.geoFrame.n, wo);
					//if (wiDotGeoN * Frame::cosTheta( its.toLocal(wi) ) <= 0 ||
					//	woDotGeoN * Frame::cosTheta( its.toLocal(wo) ) <= 0) {
					//	break;                
     //               }

					/* Keep track of the weight, medium and relative
					   refractive index along the path */
					throughput *= bsdfWeight;
                    eta *= gsampler.getEta();
					if ( its.isMediumTransition() ) {
						medium = its.getTargetMedium(woDotGeoN);
                    }

					if (gsampler.getLastSampledComponent() & BSDF::ENull) {
						++nullInteractions;
                    } else {
						delta = gsampler.getLastSampledComponent() & BSDF::EDelta;
                    }

					ray.setOrigin(its.p);
					ray.setDirection(wo);
					ray.mint = Epsilon;
				}

				if (depth++ >= m_rrDepth) {
					/* Russian roulette: try to keep path weights equal to one,
					   Stop with at least some probability to avoid
					   getting stuck (e.g. due to total internal reflection) */
					//Float q = std::min(/*throughput.average()*/std::max(bsdfWeight.max(),albedo.max()), (Float) /*0.95f*/ 1.f );
                    Float q = russianRoulette( throughput, power, eta, 
                                               m_ww, m_cfg.m_mitsuba.useWeightWindow );
					if (m_sampler->next1D() >= q) {
                        m_workResult->m_particleStats.incKills();
						break;                    
                    }
					throughput /= q;                    
				} 
			}
		}

		m_workResult->nextParticle();
		m_workResult = NULL;
	}

	void handleNewParticle() {
		m_workResult->nextParticle();
	}

	void handleSurfaceInteraction(int depth_, int nullInteractions, bool delta,
			const Intersection &its, const Medium *medium,
			const Spectrum &weight, const Float & distance) {        
        m_workResult->m_particleStats.observe( weight, m_ww.getWeightScale(), depth_ );
		int bsdfType = its.getBSDF()->getType(), depth = depth_ - nullInteractions;
		if (!(bsdfType & BSDF::EDiffuseReflection) && !(bsdfType & BSDF::EGlossyReflection))
			return;

		if ((m_type == GuidedGatherPhotonProcess::ECausticPhotons && depth > 1 && delta)
		 || (m_type == GuidedGatherPhotonProcess::ESurfacePhotons && depth > 1 && !delta)
		 || (m_type == GuidedGatherPhotonProcess::EAllIndirectSurfacePhotons && depth > 1)
		 || (m_type == GuidedGatherPhotonProcess::EAllSurfacePhotons)) {
			Photon photon(its.p, its.geoFrame.n, -its.toWorld(its.wi), weight, depth);
			photon.setDistance( distance );
			m_workResult->put( photon );
		}
	}

	void handleBackgroundInteraction( int depth, const Ray & ray, const Spectrum & weight ) {
		if ( depth > 1 ) {
			m_workResult->putBgParticle( Photon( ray.o, Normal(), ray.d, weight, depth ) );
		}
	}

	void handleMediumInteraction(int depth, int nullInteractions, bool delta,
			const MediumSamplingRecord &mRec, const Medium *medium,
			const Vector &wi, const Spectrum &weight) {
		if (m_type == GuidedGatherPhotonProcess::EVolumePhotons)
			m_workResult->put(Photon(mRec.p, Normal(0.0f, 0.0f, 0.0f),
				-wi, weight, depth-nullInteractions));
	}

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~GuidedGatherPhotonWorker() { }
protected:
	GuidedGatherPhotonProcess::EGatherType m_type;
	GuidedGatherPhotonProcess::ETraceDirection m_traceDirection;
	Importance::Sampler * m_importance;	
	size_t m_granularity;
	ref<GuidedPhotonVector> m_workResult;
	int m_maxPathDepth;
    const WeightWindow m_ww;
    const GuidingConfig m_cfg;
};

GuidedGatherPhotonProcess::GuidedGatherPhotonProcess(EGatherType type, size_t photonCount,
	size_t granularity, int maxDepth, int rrDepth, 
    const WeightWindow & ww,
    bool isLocal, bool autoCancel,
	const void *progressReporterPayload,
	Importance::Sampler * importance,	
	ETraceDirection traceDirection,
    const GuidingConfig & cfg,
	size_t numGenerated
	)
	: ParticleProcess(ParticleProcess::ETrace, photonCount, granularity, traceDirection == EFromEmitters ? "Gathering photons" : "Gathering importons",
	  progressReporterPayload, numGenerated), m_type(type), m_photonCount(photonCount), m_maxDepth(maxDepth),
	  m_rrDepth(rrDepth),  m_isLocal(isLocal), m_autoCancel(autoCancel), m_excess(0), m_numShot(0),
	  m_importance(importance),
	  m_traceDirection(traceDirection),
      m_ww( ww ),
      m_cfg( cfg ) {
	m_photonMap = new PhotonMap( photonCount * std::min( 40, std::max( 1, maxDepth ) ) );
	m_bgParticles   = new BgParticlesVec();
}

bool GuidedGatherPhotonProcess::isLocal() const {
	return m_isLocal;
}

ref<WorkProcessor> GuidedGatherPhotonProcess::createWorkProcessor() const {
	return new GuidedGatherPhotonWorker(m_type, m_granularity, m_maxDepth,
		m_rrDepth, m_ww, m_importance, m_cfg, m_traceDirection);
}

void GuidedGatherPhotonProcess::processResult(const WorkResult *wr, bool cancelled) {
	if (cancelled)
		return;
	const GuidedPhotonVector &vec = *static_cast<const GuidedPhotonVector *>(wr);
	LockGuard lock(m_resultMutex);

	/* process background particles */    
	m_bgParticles->m_data.insert( m_bgParticles->m_data.begin(), vec.bgParticlesBegin(), vec.bgParticlesEnd() );

	size_t nParticles = 0;
	for (size_t i=0; i<vec.getParticleCount(); ++i) {
		size_t start = vec.getParticleIndex(i),
			   end   = vec.getParticleIndex(i+1);
		++nParticles;
		bool full = false;
		for (size_t j=start; j<end; ++j) {
			if (!m_photonMap->tryAppend(vec[j])) {
				m_excess += vec.size() - j;
				full = true;
				break;
			}
		}
		if (full)
			break;
	}
	m_numShot += nParticles;
    m_particleStats += vec.m_particleStats;
	increaseResultCount(vec.size());
}

ParallelProcess::EStatus GuidedGatherPhotonProcess::generateWork(WorkUnit *unit, int worker) {
	/* Use the same approach as PBRT for auto canceling */
	//LockGuard lock(m_resultMutex);
 //   if (m_autoCancel && m_numShot > 100000
 //       && unsuccessful(m_photonCount, m_photonMap->size(), m_numShot)) {
 //           Log(EInfo, "Not enough photons could be collected, giving up");
 //           return EFailure;
 //   }

	return ParticleProcess::generateWork(unit, worker);
}

MTS_IMPLEMENT_CLASS(GuidedGatherPhotonProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS_S(GuidedGatherPhotonWorker, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(GuidedPhotonVector, false, WorkResult)
MTS_NAMESPACE_END
