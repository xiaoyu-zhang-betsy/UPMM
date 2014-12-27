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


#include <mitsuba/render/guided_particletracing.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/range.h>
#include <mitsuba/render/guided_brdf.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief This work result implementation stores a sequence of particles, which can be
 * sent over the wire as needed.
 *
 * It is used to implement parallel networked particle tracing passes.
 */
class GuidedParticleVector : public WorkResult {
public:
	GuidedParticleVector() { }

	inline void nextParticle() {
		m_particleIndices.push_back((uint32_t) m_particles.size());
	}

	inline void put(const Photon &p) {
		m_particles.push_back(p);
	}

	inline size_t size() const {
		return m_particles.size();
	}

	inline size_t getParticleCount() const {
		return m_particleIndices.size()-1;
	}

	inline size_t getParticleIndex(size_t idx) const {
		return m_particleIndices.at(idx);
	}

	inline void clear() {
		m_particles.clear();
		m_particleIndices.clear();
		m_bgParticles.m_data.clear();
	}

	inline const Photon &operator[](size_t index) const {
		return m_particles[index];
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
		oss << "GuidedParticleVector[size=" << m_particles.size() << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
public:
    ParticleStats m_particleStats;
protected:
	// Virtual destructor
	virtual ~GuidedParticleVector() { }
private:
	std::vector<Photon> m_particles;
	std::vector<uint32_t> m_particleIndices;
	BgParticlesVec m_bgParticles;    
};

/**
 * This class does the actual particle tracing work
 */
class GuidedEmitParticleWorker : public ParticleTracer {
protected:
    struct PathContext {
        Ray ray;
        Intersection its;
        int depth;
        const Medium * medium;
        Spectrum throughput;
        Spectrum power;
        int nullInteractions;
        bool delta;
        Float eta;
        MediumSamplingRecord mRec;
    };

    typedef std::stack<PathContext,std::vector<PathContext>> PathContextStack;
public:
	GuidedEmitParticleWorker(
        GuidedEmitParticleProcess::EGatherType type, 
        size_t granularity,
		int maxDepth, 
        int rrDepth, 
        const WeightWindow & ww,
        Importance::Sampler * importance,		
        const GuidingConfig & cfg,
        GuidedEmitParticleProcess::ETraceDirection traceDirection) 
        : ParticleTracer(maxDepth, rrDepth, false),
		m_type(type), 
		m_traceDirection( traceDirection ), 
		m_importance( importance ), 
		m_granularity(granularity), 		
        m_ww( ww ),
        m_cfg( cfg ) {}

	GuidedEmitParticleWorker(Stream *stream, InstanceManager *manager)
	 : ParticleTracer(stream, manager) {
		throw std::runtime_error("Serialization is not supported!");
	}

	ref<WorkProcessor> clone() const {
		return new GuidedEmitParticleWorker(m_type, m_granularity, m_maxDepth,
			m_rrDepth, m_ww, m_importance, m_cfg, m_traceDirection);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		throw std::runtime_error("Serialization is not supported!");
	}

	ref<WorkResult> createWorkResult() const {
		return new GuidedParticleVector();
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult,
		const bool &stop) {
            /* Warn: do not use this process() implementation that contains possible splitting with QMC samplers */

		m_workResult = static_cast<GuidedParticleVector *>(workResult);
		m_workResult->clear();

		const RangeWorkUnit *range = static_cast<const RangeWorkUnit *>(workUnit);		

		ref<Sensor> sensor    = m_scene->getSensor();
		bool needsTimeSample  = sensor->needsTimeSample();
		PositionSamplingRecord pRec(sensor->getShutterOpen()
			+ 0.5f * sensor->getShutterOpenTime());		

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
            PathContext pc;
			Point2 pixel;

			if ( m_traceDirection == GuidedEmitParticleProcess::EFromEmitters ) {                
				/* Sample both components together, which is potentially
					faster / uses a better sampling strategy */
				pc.power = m_scene->sampleEmitterRay(pc.ray, emitter,
					m_sampler->next2D(), m_sampler->next2D(), pRec.time);
				pc.medium = emitter->getMedium();
                if ( pc.power.isZero() ) {
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
				pc.power = sensor->sampleRay( pc.ray, pixel, apertureSample, timeSample );
				pc.medium = sensor->getMedium();
				handleNewParticle();
			}

			pc.depth               = 1;
            pc.nullInteractions    = 0;
			pc.delta               = false;
            pc.eta                 = 1.f;
            pc.throughput          = Spectrum( 1.0f ); // unitless path throughput (used for russian roulette)

            PathContextStack pstack;                 
            handleSplitting( m_scene->rayIntersect( pc.ray, pc.its ), pc, pstack );
            pstack.push( pc );
            
            while( !pstack.empty() ) {			
                const PathContext top = pstack.top();
                pstack.pop();
                traceParticle( top, pstack );
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

		if ((m_type == GuidedEmitParticleProcess::ECausticParticles && depth > 1 && delta)
		 || (m_type == GuidedEmitParticleProcess::ESurfaceParticles && depth > 1 && !delta)
		 || (m_type == GuidedEmitParticleProcess::EAllIndirectSurfaceParticles && depth > 1)
		 || (m_type == GuidedEmitParticleProcess::EAllSurfaceParticles)) {
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
		if (m_type == GuidedEmitParticleProcess::EVolumeParticles)
			m_workResult->put(Photon(mRec.p, Normal(0.0f, 0.0f, 0.0f),
				-wi, weight, depth-nullInteractions));
	}

    bool isPurelyDelta( const BSDF * bsdf ) const {
        return ( bsdf->getType() & BSDF::EAll & ~BSDF::EDelta ) == 0;
    }

	MTS_DECLARE_CLASS()
protected:    
	/// Virtual destructor
	virtual ~GuidedEmitParticleWorker() { }

    inline void handleSplitting( bool intersected, PathContext & pc, PathContextStack & pstack ) {        
        /* Splitting: do not split if there is no non-delta intersection in the sampled direction 
            and we are not in a medium */
        if ( m_cfg.m_mitsuba.useWeightWindow && 
            ( ( intersected && !isPurelyDelta( pc.its.getBSDF() ) ) || pc.medium ) && 
            (pc.depth < m_maxDepth || m_maxDepth < 0) && 
            m_ww.shouldSplit( pc.throughput * pc.power ) ) {

            const int count = m_ww.split( pc.throughput * pc.power );
            pc.throughput /= (Float) count;
            for ( int c = 1; c < count; c++ ) {
                pstack.push( pc );
            }

            m_workResult->m_particleStats.incSplit();
            m_workResult->m_particleStats.incNewBySplit( count );
        }
    }

    inline void traceParticle( PathContext pc, PathContextStack & pathStack ) {                
		while (!pc.throughput.isZero() && (pc.depth <= m_maxDepth || m_maxDepth < 0)) {			
			Spectrum albedo;
			Spectrum bsdfWeight;                        

			/* ==================================================================== */
			/*                 Radiative Transfer Equation sampling                 */
			/* ==================================================================== */            
			if (pc.medium && pc.medium->sampleDistance(Ray(pc.ray, 0, pc.its.t), pc.mRec, m_sampler)) {
				/* Sample the integral
					\int_x^y tau(x, x') [ \sigma_s \int_{S^2} \rho(\omega,\omega') L(x,\omega') d\omega' ] dx'
				*/

				pc.throughput *= pc.mRec.sigmaS * pc.mRec.transmittance / pc.mRec.pdfSuccess;

				/* Forward the medium scattering event to the attached handler */
				handleMediumInteraction(pc.depth, pc.nullInteractions,
						pc.delta, pc.mRec, pc.medium, -pc.ray.d, pc.throughput*pc.power);

				PhaseFunctionSamplingRecord pRec(pc.mRec, -pc.ray.d, EImportance);

				pc.throughput *= pc.medium->getPhaseFunction()->sample(pRec, m_sampler);
				pc.delta = false;

				pc.ray = Ray(pc.mRec.p, pRec.wo, pc.ray.time);
				pc.ray.mint = 0;
			} else if (pc.its.t == std::numeric_limits<Float>::infinity()) {
				/* There is no surface in this direction */
				if ( m_scene->hasEnvironmentEmitter() && m_traceDirection == GuidedEmitParticleProcess::EFromSensor ) {
					if ( pc.medium ) {
						pc.throughput *= pc.mRec.transmittance / pc.mRec.pdfFailure;
                    }
					handleBackgroundInteraction( pc.depth, pc.ray, pc.throughput * pc.power );
				}
				break;
			} else {
				/* Sample
					tau(x, y) (Surface integral). This happens with probability mRec.pdfFailure
					Account for this and multiply by the proper per-color-channel transmittance.
				*/
				if ( pc.medium ) {
					pc.throughput *= pc.mRec.transmittance / pc.mRec.pdfFailure;
                }
					
				Float distance = (pc.its.p - pc.ray.o).length();
					
                if ( pc.throughput.max() < 9e-7f ) {
                    __debugbreak();
                }

				/* Forward the surface scattering event to the attached handler */
				handleSurfaceInteraction(pc.depth, pc.nullInteractions, pc.delta, 
                    pc.its, pc.medium, pc.throughput * pc.power, distance);

				/* Prepare importance & bsdf sampling distribution */
				RayDifferential rayd( pc.ray );
				GuidedBRDF gsampler( pc.its, rayd, m_importance, m_cfg.m_mitsuba.bsdfSamplingProbability );

				/* Sample ray direction */
				Vector wo;
				Float pdf;
				Spectrum albedo;
				bsdfWeight = gsampler.sample(wo, pdf, m_sampler, albedo);

				if ( bsdfWeight.isZero() || pdf == 0.0f ) {
					break;
                }

				/* Prevent light leaks due to the use of shading normals -- [Veach, p. 158] */
				Vector wi = -pc.ray.d;
				//Float wiDotGeoN = dot(pc.its.geoFrame.n, wi);
				Float woDotGeoN = dot(pc.its.geoFrame.n, wo);
				//if (wiDotGeoN * Frame::cosTheta( pc.its.toLocal(wi) ) <= 0 ||
				//	woDotGeoN * Frame::cosTheta( pc.its.toLocal(wo) ) <= 0) {
				//	break;                
    //               }                

				/* Keep track of the weight, medium and relative
					refractive index along the path */
				pc.throughput *= bsdfWeight;
                
                pc.eta *= gsampler.getEta();
				if ( pc.its.isMediumTransition() ) {
					pc.medium = pc.its.getTargetMedium( woDotGeoN );
                }

				if (gsampler.getLastSampledComponent() & BSDF::ENull) {
					++pc.nullInteractions;
                } else {
					pc.delta = gsampler.getLastSampledComponent() & BSDF::EDelta;
                }

				pc.ray.setOrigin( pc.its.p );
				pc.ray.setDirection( wo );
				pc.ray.mint = Epsilon;
			}

			if (pc.depth++ >= m_rrDepth) {
				/* Russian roulette: try to keep path weights equal to one,
					Stop with at least some probability to avoid
					getting stuck (e.g. due to total internal reflection) */
                Float q = russianRoulette( pc.throughput, pc.power, pc.eta, 
                                            m_ww, m_cfg.m_mitsuba.useWeightWindow );
				if ( m_sampler->next1D() >= q ) {
                    m_workResult->m_particleStats.incKills();
					break;                    
                }
				pc.throughput /= q;                    
			}            
                
            handleSplitting( m_scene->rayIntersect( pc.ray, pc.its ), pc, pathStack );
		}
    }
protected:
	GuidedEmitParticleProcess::EGatherType m_type;
	GuidedEmitParticleProcess::ETraceDirection m_traceDirection;
	Importance::Sampler * m_importance;	
	size_t m_granularity;
	ref<GuidedParticleVector> m_workResult;
	int m_maxPathDepth;
    const WeightWindow m_ww;
    const GuidingConfig m_cfg;
};

GuidedEmitParticleProcess::GuidedEmitParticleProcess(EGatherType type, size_t photonCount,
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

bool GuidedEmitParticleProcess::isLocal() const {
	return m_isLocal;
}

ref<WorkProcessor> GuidedEmitParticleProcess::createWorkProcessor() const {
	return new GuidedEmitParticleWorker(m_type, m_granularity, m_maxDepth,
		m_rrDepth, m_ww, m_importance, m_cfg, m_traceDirection);
}

void GuidedEmitParticleProcess::processResult(const WorkResult *wr, bool cancelled) {
	if (cancelled)
		return;
	const GuidedParticleVector &vec = *static_cast<const GuidedParticleVector *>(wr);
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

ParallelProcess::EStatus GuidedEmitParticleProcess::generateWork(WorkUnit *unit, int worker) {
	/* Use the same approach as PBRT for auto canceling */
	//LockGuard lock(m_resultMutex);
 //   if (m_autoCancel && m_numShot > 100000
 //       && unsuccessful(m_photonCount, m_photonMap->size(), m_numShot)) {
 //           Log(EInfo, "Not enough photons could be collected, giving up");
 //           return EFailure;
 //   }

	return ParticleProcess::generateWork(unit, worker);
}

MTS_IMPLEMENT_CLASS(GuidedEmitParticleProcess, false, ParticleProcess)
MTS_IMPLEMENT_CLASS_S(GuidedEmitParticleWorker, false, ParticleTracer)
MTS_IMPLEMENT_CLASS(GuidedParticleVector, false, WorkResult)
MTS_NAMESPACE_END
