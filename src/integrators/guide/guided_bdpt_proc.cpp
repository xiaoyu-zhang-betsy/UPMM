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

#include <mitsuba/render/guiding.h>
#include <mitsuba/render/guided_brdf.h>

#include "guided_bdpt_proc.h"
#include <mitsuba/core/statistics.h>
#include <mitsuba/core/sfcurve.h>
#include <mitsuba/bidir/util.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Worker implementation                        */
/* ==================================================================== */

class GuidedBDPTRenderer : public WorkProcessor {
public:
	GuidedBDPTRenderer(const GuidedBDPTConfiguration &config) : m_config(config) { }

	GuidedBDPTRenderer(Stream *stream, InstanceManager *manager)
		: WorkProcessor(stream, manager), m_config(stream) { }

	virtual ~GuidedBDPTRenderer() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		m_config.serialize(stream);
	}

	ref<WorkUnit> createWorkUnit() const {
		return new RectangularWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new GuidedBDPTWorkResult(m_config, m_rfilter.get(),
			Vector2i(m_config.blockSize));
	}

	void prepare() {
		Scene *scene = static_cast<Scene *>(getResource("scene"));
		m_scene = new Scene(scene);
		m_sampler = static_cast<Sampler *>(getResource("sampler"));
		m_sensor = static_cast<Sensor *>(getResource("sensor"));
		m_rfilter = m_sensor->getFilm()->getReconstructionFilter();
		m_scene->removeSensor(scene->getSensor());
		m_scene->addSensor(m_sensor);
		m_scene->setSensor(m_sensor);
		m_scene->setSampler(m_sampler);
		m_scene->wakeup(NULL, m_resources);
		m_scene->initializeBidirectional();

		m_guidingSampler = static_cast<GuidingSamplers *>(getResource("guidingSampler"));
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
		const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(workUnit);
		GuidedBDPTWorkResult *result = static_cast<GuidedBDPTWorkResult *>(workResult);
		bool needsTimeSample = m_sensor->needsTimeSample();
		Float time = m_sensor->getShutterOpen();

		result->setOffset(rect->getOffset());
		result->setSize(rect->getSize());
		result->clear();
		m_hilbertCurve.initialize(TVector2<uint8_t>(rect->getSize()));

		#if defined(MTS_DEBUG_FP)
			enableFPExceptions();
		#endif

		Path emitterSubpath;
		Path sensorSubpath;

		/* Determine the necessary random walk depths based on properties of
		   the endpoints */
		int emitterDepth = m_config.maxDepth,
		    sensorDepth = m_config.maxDepth;

		/* Go one extra step if the sensor can be intersected */
		if (!m_scene->hasDegenerateSensor() && emitterDepth != -1)
			++emitterDepth;

		/* Go one extra step if there are emitters that can be intersected */
		if (!m_scene->hasDegenerateEmitters() && sensorDepth != -1)
			++sensorDepth;

		for (size_t i=0; i<m_hilbertCurve.getPointCount(); ++i) {
			Point2i offset = Point2i(m_hilbertCurve[i]) + Vector2i(rect->getOffset());
			m_sampler->generate(offset);

			for (size_t j = 0; j<m_sampler->getSampleCount(); j++) {
				if (stop)
					break;

				if (needsTimeSample)
					time = m_sensor->sampleTime(m_sampler->next1D());

				/* Start new emitter and sensor subpaths */
				emitterSubpath.initialize(m_scene, time, EImportance, m_pool);
				sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);

				/* Perform a random walk using alternating steps on each path */
				alternatingRandomWalkFromPixel(m_scene, m_sampler,
					emitterSubpath, emitterDepth, sensorSubpath,
					sensorDepth, offset, m_config.rrDepth, m_pool);

				evaluate(result, emitterSubpath, sensorSubpath);

				emitterSubpath.release(m_pool);
				sensorSubpath.release(m_pool);

				m_sampler->advance();
			}
		}

		#if defined(MTS_DEBUG_FP)
			disableFPExceptions();
		#endif

		/* Make sure that there were no memory leaks */
		Assert(m_pool.unused());
	}

	bool sampleNext(PathVertex* current,
		const Scene *scene, Sampler *sampler,
		const PathVertex *pred, const PathEdge *predEdge,
		PathEdge *succEdge, PathVertex *succ,
		ETransportMode mode, bool russianRoulette, Spectrum *throughput) {
		Ray ray;

		memset(succEdge, 0, sizeof(PathEdge));
		memset(succ, 0, sizeof(PathVertex));

		succEdge->medium = (predEdge == NULL) ? NULL : predEdge->medium;
		current->rrWeight = 1.0f;

		switch (current->type) {
		case PathVertex::EEmitterSupernode: {
			BDAssert(mode == EImportance && pred == NULL && predEdge == NULL);
			PositionSamplingRecord &pRec = succ->getPositionSamplingRecord();
			const EndpointRecord &eRec = current->getEndpointRecord();
			pRec = PositionSamplingRecord(eRec.time);
			Spectrum result = scene->sampleEmitterPosition(pRec, sampler->next2D());
			if (result.isZero())
				return false;

			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			current->weight[EImportance] = result;
			current->pdf[EImportance] = pRec.pdf;
			current->measure = pRec.measure;
			succ->type = current->EEmitterSample;
			succ->degenerate = emitter->getType() & Emitter::EDeltaDirection;

			succEdge->weight[EImportance] = Spectrum(1.0f);
			succEdge->pdf[EImportance] = 1.0f;
			succEdge->medium = emitter->getMedium();

			return true;
		}
			break;

		case PathVertex::ESensorSupernode: {
			BDAssert(mode == ERadiance && pred == NULL && predEdge == NULL);
			PositionSamplingRecord &pRec = succ->getPositionSamplingRecord();
			const EndpointRecord &eRec = current->getEndpointRecord();
			pRec = PositionSamplingRecord(eRec.time);
			Spectrum result = scene->sampleSensorPosition(pRec, sampler->next2D());
			if (result.isZero())
				return false;

			const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
			current->weight[ERadiance] = result;
			current->pdf[ERadiance] = pRec.pdf;
			current->measure = pRec.measure;
			succ->type = PathVertex::ESensorSample;
			succ->degenerate = sensor->getType()
				& Sensor::EDeltaDirection;

			succEdge->weight[ERadiance] = Spectrum(1.0f);
			succEdge->pdf[ERadiance] = 1.0f;
			succEdge->medium = sensor->getMedium();

			return true;
		}
			break;

		case PathVertex::EEmitterSample: {
			BDAssert(mode == EImportance && pred->type == PathVertex::EEmitterSupernode);
			PositionSamplingRecord &pRec = current->getPositionSamplingRecord();
			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			DirectionSamplingRecord dRec;

			Spectrum result = emitter->sampleDirection(dRec, pRec,
				emitter->needsDirectionSample() ? sampler->next2D() : Point2(0.5f));

			if (result.isZero())
				return false;

			current->weight[EImportance] = result;
			current->weight[ERadiance] = result * dRec.pdf * (
				emitter->isOnSurface() ? 1.0f / absDot(dRec.d, pRec.n) : 1.0f);
			current->pdf[EImportance] = dRec.pdf;
			current->pdf[ERadiance] = 1.0f;

			current->measure = dRec.measure;
			succEdge->medium = emitter->getMedium();
			ray.time = pRec.time;
			ray.setOrigin(pRec.p);
			ray.setDirection(dRec.d);
		}
			break;

		case PathVertex::ESensorSample: {
			BDAssert(mode == ERadiance && pred->type == PathVertex::ESensorSupernode);
			PositionSamplingRecord &pRec = current->getPositionSamplingRecord();
			const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
			DirectionSamplingRecord dRec;

			Spectrum result = sensor->sampleDirection(dRec, pRec,
				sensor->needsDirectionSample() ? sampler->next2D() : Point2(0.5f));

			if (result.isZero())
				return false;

			current->weight[EImportance] = result * dRec.pdf * (
				sensor->isOnSurface() ? 1.0f / absDot(dRec.d, pRec.n) : 1.0f);
			current->weight[ERadiance] = result;
			current->pdf[EImportance] = 1.0f;
			current->pdf[ERadiance] = dRec.pdf;

			current->measure = dRec.measure;
			succEdge->medium = sensor->getMedium();
			ray.time = pRec.time;
			ray.setOrigin(pRec.p);
			ray.setDirection(dRec.d);
		}
			break;

		case PathVertex::ESurfaceInteraction: {
			Intersection its = current->getIntersection();
			//const BSDF *bsdf = its.getBSDF();
			Vector wi = normalize(pred->getPosition() - its.p);
			Vector wo;

			GuidedBRDF gsampler(its, (mode == ERadiance) ? m_guidingSampler->getRadianceSampler() : m_guidingSampler->getImportanceSampler(),
				m_guidingSampler->getConfig().m_mitsuba.bsdfSamplingProbability);					
			if (!gsampler.isValid()) {
				m_guidingSampler->distributionConstructionFailed(its);
			}

			/* Sample the BSDF */
			BSDFSamplingRecord bRec(its, sampler, mode);
			bRec.wi = its.toLocal(wi);
			//current->weight[mode] = bsdf->sample(bRec, current->pdf[mode], sampler->next2D());
			current->weight[mode] = gsampler.sample(wo, current->pdf[mode], sampler);
			if (current->weight[mode].isZero())
				return false;

			bRec.sampledType = gsampler.getLastSampledComponent();
			bRec.wo = its.toLocal(wo);
			bRec.eta = gsampler.getEta();
			if (bRec.sampledType == 0) bRec.sampledType |= BSDF::ESmooth;


			current->measure = BSDF::getMeasure(bRec.sampledType);
			current->componentType = (uint16_t)(bRec.sampledType & BSDF::EAll);

			wo = its.toWorld(bRec.wo);

			/* Prevent light leaks due to the use of shading normals */
			Float wiDotGeoN = dot(its.geoFrame.n, wi),
				woDotGeoN = dot(its.geoFrame.n, wo);
			if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
				woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				return false;

			/* Compute the reverse quantities */
			bRec.reverse();
			Intersection its_inv = its;
			its_inv.wi = bRec.wi;
			GuidedBRDF gsampler_inv(its_inv, (mode == ERadiance) ? m_guidingSampler->getImportanceSampler() : m_guidingSampler->getRadianceSampler(),
				m_guidingSampler->getConfig().m_mitsuba.bsdfSamplingProbability);
			current->pdf[1 - mode] = gsampler_inv.pdf(its.toWorld(bRec.wo));
			//current->pdf[1 - mode] = bsdf->pdf(bRec, (EMeasure)(current->measure));
			if (current->pdf[1 - mode] == 0) {
				/* This can happen rarely due to roundoff errors -- be strict */
				return false;
			}
			
			//if (!(bsdf->getType() & BSDF::ENonSymmetric)) {
			if (!gsampler.isBSDFPureSpecular() && false) {
				/* Make use of symmetry -- no need to re-evaluate
				everything (only the pdf and cosine factors changed) */
				current->weight[1 - mode] = current->weight[mode] * (current->pdf[mode] / current->pdf[1 - mode]);
				if (current->measure == ESolidAngle)
					current->weight[1 - mode] *=
					std::abs(Frame::cosTheta(bRec.wo) / Frame::cosTheta(bRec.wi));
			}
			else {
				current->weight[1 - mode] = gsampler_inv.eval(its.toWorld(bRec.wo)) / current->pdf[1 - mode];  
				//current->weight[1 - mode] = bsdf->eval(bRec, (EMeasure)(current->measure)) / current->pdf[1 - mode];
			}
			bRec.reverse();

			/* Adjoint BSDF for shading normals */
			if (mode == EImportance)
				current->weight[EImportance] *= std::abs(
				(Frame::cosTheta(bRec.wi) * woDotGeoN) /
				(Frame::cosTheta(bRec.wo) * wiDotGeoN));
			else
				current->weight[EImportance] *= std::abs(
				(Frame::cosTheta(bRec.wo) * wiDotGeoN) /
				(Frame::cosTheta(bRec.wi) * woDotGeoN));

			/// For BDPT & russian roulette, track radiance * eta^2
			if (throughput && mode == ERadiance && bRec.eta != 1)
				(*throughput) *= bRec.eta * bRec.eta;

			ray.time = its.time;
			ray.setOrigin(its.p);
			ray.setDirection(wo);
		}
			break;

		default:
			SLog(EError, "PathVertex::sampleNext(): Encountered an "
				"unsupported vertex type (%i)!", current->type);
			return false;
		}

		if (throughput) {
			/* For BDPT: keep track of the path throughput to run russian roulette */
			(*throughput) *= current->weight[mode];

			if (russianRoulette) {
				Float q = std::min(throughput->max(), (Float) 0.95f);

				if (sampler->next1D() > q) {
					current->measure = EInvalidMeasure;
					return false;
				}
				else {
					current->rrWeight = 1.0f / q;
					(*throughput) *= current->rrWeight;
				}
			}
		}

		if (!succEdge->sampleNext(scene, sampler, current, ray, succ, mode)) {
			/* Sampling a successor edge + vertex failed, hence the vertex
			is not committed to a particular measure yet -- revert. */
			current->measure = EInvalidMeasure;
			return false;
		}
		else {
			if (throughput)
				(*throughput) *= succEdge->weight[mode];
		}

		/* Convert from solid angle to area measure */
		if (current->measure == ESolidAngle) {
			current->measure = EArea;

			current->pdf[mode] /= succEdge->length * succEdge->length;
			if (succ->isOnSurface())
				current->pdf[mode] *= absDot(ray.d, succ->getGeometricNormal());

			if (predEdge->length != 0.0f) {
				current->pdf[1 - mode] /= predEdge->length * predEdge->length;
				if (pred->isOnSurface())
					current->pdf[1 - mode] *= absDot(predEdge->d, pred->getGeometricNormal());
			}
		}

		return true;
	}

	void alternatingRandomWalkFromPixel(const Scene *scene, Sampler *sampler,
		Path &emitterPath, int nEmitterSteps, Path &sensorPath, int nSensorSteps,
		const Point2i &pixelPosition, int rrStart, MemoryPool &pool) {
		/* Determine the relevant edges and vertices to start the random walk */
		PathVertex *curVertexS = emitterPath.vertex(0),
			*curVertexT = sensorPath.vertex(0),
			*predVertexS = NULL, *predVertexT = NULL;
		PathEdge   *predEdgeS = NULL, *predEdgeT = NULL;

		PathVertex *v1 = pool.allocVertex(), *v2 = pool.allocVertex();
		PathEdge *e0 = pool.allocEdge(), *e1 = pool.allocEdge();

		/* Use a special sampling routine for the first two sensor vertices so that
		the resulting subpath passes through the specified pixel position */
		int t = curVertexT->sampleSensor(scene,
			sampler, pixelPosition, e0, v1, e1, v2);

		if (t >= 1) {
			sensorPath.append(e0, v1);
		}
		else {
			pool.release(e0);
			pool.release(v1);
		}

		if (t == 2) {
			sensorPath.append(e1, v2);
			predVertexT = v1;
			curVertexT = v2;
			predEdgeT = e1;
		}
		else {
			pool.release(e1);
			pool.release(v2);
			curVertexT = NULL;
		}

		Spectrum throughputS(1.0f), throughputT(1.0f);

		int s = 0;
		do {
			if (curVertexT && (t < nSensorSteps || nSensorSteps == -1)) {
				PathVertex *succVertexT = pool.allocVertex();
				PathEdge *succEdgeT = pool.allocEdge();

				if (sampleNext(curVertexT,
					scene, sampler, predVertexT,
					predEdgeT, succEdgeT, succVertexT, ERadiance,
					rrStart != -1 && t >= rrStart, &throughputT)) {
					sensorPath.append(succEdgeT, succVertexT);
					predVertexT = curVertexT;
					curVertexT = succVertexT;
					predEdgeT = succEdgeT;
					t++;
				}
				else {
					pool.release(succVertexT);
					pool.release(succEdgeT);
					curVertexT = NULL;
				}

// 				if (curVertexT->sampleNext(scene, sampler, predVertexT,
// 					predEdgeT, succEdgeT, succVertexT, ERadiance,
// 					rrStart != -1 && t >= rrStart, &throughputT)) {
// 					sensorPath.append(succEdgeT, succVertexT);
// 					predVertexT = curVertexT;
// 					curVertexT = succVertexT;
// 					predEdgeT = succEdgeT;
// 					t++;
// 				}
// 				else {
// 					pool.release(succVertexT);
// 					pool.release(succEdgeT);
// 					curVertexT = NULL;
// 				}

			}
			else {
				curVertexT = NULL;
			}

			if (curVertexS && (s < nEmitterSteps || nEmitterSteps == -1)) {
				PathVertex *succVertexS = pool.allocVertex();
				PathEdge *succEdgeS = pool.allocEdge();

				if (sampleNext(curVertexS,
					scene, sampler, predVertexS,
					predEdgeS, succEdgeS, succVertexS, EImportance,
					rrStart != -1 && s >= rrStart, &throughputS)) {
					emitterPath.append(succEdgeS, succVertexS);
					predVertexS = curVertexS;
					curVertexS = succVertexS;
					predEdgeS = succEdgeS;
					s++;
				}
				else {
					pool.release(succVertexS);
					pool.release(succEdgeS);
					curVertexS = NULL;
				}

// 				if (curVertexS->sampleNext(scene, sampler, predVertexS,
// 					predEdgeS, succEdgeS, succVertexS, EImportance,
// 					rrStart != -1 && s >= rrStart, &throughputS)) {
// 					emitterPath.append(succEdgeS, succVertexS);
// 					predVertexS = curVertexS;
// 					curVertexS = succVertexS;
// 					predEdgeS = succEdgeS;
// 					s++;
// 				}
// 				else {
// 					pool.release(succVertexS);
// 					pool.release(succEdgeS);
// 					curVertexS = NULL;
// 				}
			}
			else {
				curVertexS = NULL;
			}
		} while (curVertexS || curVertexT);
	}	

	/// Evaluate the contributions of the given eye and light paths
	void evaluate(GuidedBDPTWorkResult *wr,
			Path &emitterSubpath, Path &sensorSubpath) {
		Point2 initialSamplePos = sensorSubpath.vertex(1)->getSamplePosition();
		const Scene *scene = m_scene;
		PathVertex tempEndpoint, tempSample;
		PathEdge tempEdge, connectionEdge;

		/* Compute the combined weights along the two subpaths */
		Spectrum *importanceWeights = (Spectrum *) alloca(emitterSubpath.vertexCount() * sizeof(Spectrum)),
				 *radianceWeights  = (Spectrum *) alloca(sensorSubpath.vertexCount()  * sizeof(Spectrum));

		importanceWeights[0] = radianceWeights[0] = Spectrum(1.0f);
		for (size_t i=1; i<emitterSubpath.vertexCount(); ++i)
			importanceWeights[i] = importanceWeights[i-1] *
				emitterSubpath.vertex(i-1)->weight[EImportance] *
				emitterSubpath.vertex(i-1)->rrWeight *
				emitterSubpath.edge(i-1)->weight[EImportance];

		for (size_t i=1; i<sensorSubpath.vertexCount(); ++i)
			radianceWeights[i] = radianceWeights[i-1] *
				sensorSubpath.vertex(i-1)->weight[ERadiance] *
				sensorSubpath.vertex(i-1)->rrWeight *
				sensorSubpath.edge(i-1)->weight[ERadiance];

		Spectrum sampleValue(0.0f);
		for (int s = (int) emitterSubpath.vertexCount()-1; s >= 0; --s) {
			/* Determine the range of sensor vertices to be traversed,
			   while respecting the specified maximum path length */
			int minT = std::max(2-s, m_config.lightImage ? 0 : 2),
			    maxT = (int) sensorSubpath.vertexCount() - 1;
			if (m_config.maxDepth != -1)
				maxT = std::min(maxT, m_config.maxDepth + 1 - s);

			for (int t = maxT; t >= minT; --t) {
				PathVertex
					*vsPred = emitterSubpath.vertexOrNull(s-1),
					*vtPred = sensorSubpath.vertexOrNull(t-1),
					*vs = emitterSubpath.vertex(s),
					*vt = sensorSubpath.vertex(t);
				PathEdge
					*vsEdge = emitterSubpath.edgeOrNull(s-1),
					*vtEdge = sensorSubpath.edgeOrNull(t-1);

				//if (s + t <= 3) continue; // Temporary exclude direct lighting

				RestoreMeasureHelper rmh0(vs), rmh1(vt);

				/* Will be set to true if direct sampling was used */
				bool sampleDirect = false;

				/* Stores the pixel position associated with this sample */
				Point2 samplePos = initialSamplePos;

				/* Allowed remaining number of ENull vertices that can
				   be bridged via pathConnect (negative=arbitrarily many) */
				int remaining = m_config.maxDepth - s - t + 1;

				/* Will receive the path weight of the (s, t)-connection */
				Spectrum value;

				/* Account for the terms of the measurement contribution
				   function that are coupled to the connection endpoints */
				if (vs->isEmitterSupernode()) {
					/* If possible, convert 'vt' into an emitter sample */
					if (!vt->cast(scene, PathVertex::EEmitterSample) || vt->isDegenerate())
						continue;

					value = radianceWeights[t] *
						vs->eval(scene, vsPred, vt, EImportance) *
						vt->eval(scene, vtPred, vs, ERadiance);
				} else if (vt->isSensorSupernode()) {
					/* If possible, convert 'vs' into an sensor sample */
					if (!vs->cast(scene, PathVertex::ESensorSample) || vs->isDegenerate())
						continue;

					/* Make note of the changed pixel sample position */
					if (!vs->getSamplePosition(vsPred, samplePos))
						continue;

					value = importanceWeights[s] *
						vs->eval(scene, vsPred, vt, EImportance) *
						vt->eval(scene, vtPred, vs, ERadiance);
				} else if (m_config.sampleDirect && ((t == 1 && s > 1) || (s == 1 && t > 1))) {
					/* s==1/t==1 path: use a direct sampling strategy if requested */
					if (s == 1) {
						if (vt->isDegenerate())
							continue;
						/* Generate a position on an emitter using direct sampling */
						value = radianceWeights[t] * vt->sampleDirect(scene, m_sampler,
							&tempEndpoint, &tempEdge, &tempSample, EImportance);
						if (value.isZero())
							continue;
						vs = &tempSample; vsPred = &tempEndpoint; vsEdge = &tempEdge;
						value *= vt->eval(scene, vtPred, vs, ERadiance);
						vt->measure = EArea;
					} else {
						if (vs->isDegenerate())
							continue;
						/* Generate a position on the sensor using direct sampling */
						value = importanceWeights[s] * vs->sampleDirect(scene, m_sampler,
							&tempEndpoint, &tempEdge, &tempSample, ERadiance);
						if (value.isZero())
							continue;
						vt = &tempSample; vtPred = &tempEndpoint; vtEdge = &tempEdge;
						value *= vs->eval(scene, vsPred, vt, EImportance);
						vs->measure = EArea;
					}

					sampleDirect = true;
				} else {
					/* Can't connect degenerate endpoints */
					if (vs->isDegenerate() || vt->isDegenerate())
						continue;

					value = importanceWeights[s] * radianceWeights[t] *
						vs->eval(scene, vsPred, vt, EImportance) *
						vt->eval(scene, vtPred, vs, ERadiance);

					/* Temporarily force vertex measure to EArea. Needed to
					   handle BSDFs with diffuse + specular components */
					vs->measure = vt->measure = EArea;
				}

				/* Attempt to connect the two endpoints, which could result in
				   the creation of additional vertices (index-matched boundaries etc.) */
				int interactions = remaining; // backup
				if (value.isZero() || !connectionEdge.pathConnectAndCollapse(
						scene, vsEdge, vs, vt, vtEdge, interactions))
					continue;

				/* Determine the pixel sample position when necessary */
				if (vt->isSensorSample() && !vt->getSamplePosition(vs, samplePos))
					continue;				

				/* Account for the terms of the measurement contribution
				   function that are coupled to the connection edge */
				if (!sampleDirect)
					value *= connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
				else
					value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
							(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

				if (sampleDirect) {
					/* A direct sampling strategy was used, which generated
					   two new vertices at one of the path ends. Temporarily
					   modify the path to reflect this change */
					if (t == 1)
						sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
					else
						emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
				}

				/* Compute the multiple importance sampling weight */
				Float weight = miWeight(scene, emitterSubpath, &connectionEdge,
					sensorSubpath, s, t, m_config.sampleDirect, m_config.lightImage);

				if (sampleDirect) {
					/* Now undo the previous change */
					if (t == 1)
						sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
					else
						emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
				}				
				
				#if GBDPT_DEBUG == 1
					/* When the debug mode is on, collect samples
					   separately for each sampling strategy. Note: the
					   following piece of code artificially increases the
					   exposure of longer paths */
					Spectrum splatValue = value * (m_config.showWeighted
 						? weight : 1.0f);// * std::pow(2.0f, s+t-3.0f));
					wr->putDebugSample(s, t, samplePos, splatValue);
					wr->putDebugSampleM(s, t, samplePos, value);
				#endif

					value *= weight;
				if (t >= 2)
					sampleValue += value;
				else
					wr->putLightSample(samplePos, value);
			}
		}
		wr->putSample(initialSamplePos, sampleValue);
	}	

	Float evalPdf(const PathVertex* current, const Scene *scene, const PathVertex *pred,
		const PathVertex *succ, ETransportMode mode, EMeasure measure) const {
		Vector wo(0.0f);
		Float dist = 0.0f, result = 0.0f;

		switch (current->type) {
		case PathVertex::EEmitterSupernode: {
			if (mode != EImportance || pred != NULL || succ->type != PathVertex::EEmitterSample)
				return 0.0f;
			PositionSamplingRecord pRec = succ->getPositionSamplingRecord();
			pRec.measure = measure;
			return scene->pdfEmitterPosition(pRec);
		}
			break;

		case PathVertex::ESensorSupernode: {
			if (mode != ERadiance || pred != NULL || succ->type != PathVertex::ESensorSample)
				return 0.0f;
			PositionSamplingRecord pRec = succ->getPositionSamplingRecord();
			pRec.measure = measure;
			return scene->pdfSensorPosition(pRec);
		}
			break;

		case PathVertex::EEmitterSample: {
			if (mode == ERadiance && succ->type == PathVertex::EEmitterSupernode)
				return 1.0f;
			else if (mode != EImportance || pred->type != PathVertex::EEmitterSupernode)
				return 0.0f;

			const PositionSamplingRecord &pRec = current->getPositionSamplingRecord();
			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			wo = succ->getPosition() - pRec.p;
			dist = wo.length(); wo /= dist;
			DirectionSamplingRecord dRec(wo, measure == EArea ? ESolidAngle : measure);
			result = emitter->pdfDirection(dRec, pRec);
		}
			break;

		case PathVertex::ESensorSample: {
			if (mode == EImportance && succ->type == PathVertex::ESensorSupernode)
				return 1.0f;
			else if (mode != ERadiance || pred->type != PathVertex::ESensorSupernode)
				return 0.0f;

			const PositionSamplingRecord &pRec = current->getPositionSamplingRecord();
			const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
			wo = succ->getPosition() - pRec.p;
			dist = wo.length(); wo /= dist;
			DirectionSamplingRecord dRec(wo, measure == EArea ? ESolidAngle : measure);
			result = sensor->pdfDirection(dRec, pRec);
		}
			break;

		case PathVertex::ESurfaceInteraction: {
			const Intersection &its = current->getIntersection();			
			wo = succ->getPosition() - its.p;
			dist = wo.length(); wo /= dist;			

			Point predP = pred->getPosition();
			Vector wi = normalize(predP - its.p);

			//const BSDF *bsdf = its.getBSDF();
			Intersection itsi = its;
			itsi.wi = itsi.toLocal(wi);
			GuidedBRDF gsampler(itsi, (mode == ERadiance) ? m_guidingSampler->getRadianceSampler() : m_guidingSampler->getImportanceSampler(),
				m_guidingSampler->getConfig().m_mitsuba.bsdfSamplingProbability);

			BSDFSamplingRecord bRec(its, its.toLocal(wi), its.toLocal(wo), mode);
			//result = bsdf->pdf(bRec, measure == EArea ? ESolidAngle : measure);			
			result = gsampler.pdf(wo);

			/* Prevent light leaks due to the use of shading normals */
			Float wiDotGeoN = dot(its.geoFrame.n, wi),
				woDotGeoN = dot(its.geoFrame.n, wo);

			if (wiDotGeoN * Frame::cosTheta(bRec.wi) <= 0 ||
				woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
				return 0.0f;
		}
			break;		

		default:
			SLog(EError, "PathVertex::evalPdf(): Encountered an "
				"unsupported vertex type (%i)!", current->type);
			return 0.0f;
		}

		if (measure == EArea) {
			result /= dist * dist;
			if (succ->isOnSurface())
				result *= absDot(wo, succ->getGeometricNormal());
		}

		return result;
	}

	Float miWeight(const Scene *scene, const Path &emitterSubpath,
		const PathEdge *connectionEdge, const Path &sensorSubpath,
		int s, int t, bool sampleDirect, bool lightImage) {
		int k = s + t + 1, n = k + 1;

		const PathVertex
			*vsPred = emitterSubpath.vertexOrNull(s - 1),
			*vtPred = sensorSubpath.vertexOrNull(t - 1),
			*vs = emitterSubpath.vertex(s),
			*vt = sensorSubpath.vertex(t);

		/* pdfImp[i] and pdfRad[i] store the area/volume density of vertex
		'i' when sampled from the adjacent vertex in the emitter
		and sensor direction, respectively. */

		Float ratioEmitterDirect = 0.0f, ratioSensorDirect = 0.0f;
		Float *pdfImp = (Float *)alloca(n * sizeof(Float)),
			*pdfRad = (Float *)alloca(n * sizeof(Float));
		bool  *connectable = (bool *)alloca(n * sizeof(bool)),
			*isNull = (bool *)alloca(n * sizeof(bool));

		/* Keep track of which vertices are connectable / null interactions */
		int pos = 0;
		for (int i = 0; i <= s; ++i) {
			const PathVertex *v = emitterSubpath.vertex(i);
			connectable[pos] = v->isConnectable();
			isNull[pos] = v->isNullInteraction() && !connectable[pos];
			pos++;
		}

		for (int i = t; i >= 0; --i) {
			const PathVertex *v = sensorSubpath.vertex(i);
			connectable[pos] = v->isConnectable();
			isNull[pos] = v->isNullInteraction() && !connectable[pos];
			pos++;
		}

		if (k <= 3)
			sampleDirect = false;

		EMeasure vsMeasure = EArea, vtMeasure = EArea;
		if (sampleDirect) {
			/* When direct sampling is enabled, we may be able to create certain
			connections that otherwise would have failed (e.g. to an
			orthographic camera or a directional light source) */
			const AbstractEmitter *emitter = (s > 0 ? emitterSubpath.vertex(1) : vt)->getAbstractEmitter();
			const AbstractEmitter *sensor = (t > 0 ? sensorSubpath.vertex(1) : vs)->getAbstractEmitter();

			EMeasure emitterDirectMeasure = emitter->getDirectMeasure();
			EMeasure sensorDirectMeasure = sensor->getDirectMeasure();

			connectable[0] = emitterDirectMeasure != EDiscrete && emitterDirectMeasure != EInvalidMeasure;
			connectable[1] = emitterDirectMeasure != EInvalidMeasure;
			connectable[k - 1] = sensorDirectMeasure != EInvalidMeasure;
			connectable[k] = sensorDirectMeasure != EDiscrete && sensorDirectMeasure != EInvalidMeasure;

			/* The following is needed to handle orthographic cameras &
			directional light sources together with direct sampling */
			if (t == 1)
				vtMeasure = sensor->needsDirectionSample() ? EArea : EDiscrete;
			else if (s == 1)
				vsMeasure = emitter->needsDirectionSample() ? EArea : EDiscrete;
		}

		/* Collect importance transfer area/volume densities from vertices */
		pos = 0;
		pdfImp[pos++] = 1.0;

		for (int i = 0; i<s; ++i)
			pdfImp[pos++] = emitterSubpath.vertex(i)->pdf[EImportance]
			* emitterSubpath.edge(i)->pdf[EImportance];

		pdfImp[pos++] = evalPdf(vs, scene, vsPred, vt, EImportance, vsMeasure)
			* connectionEdge->pdf[EImportance];

		if (t > 0) {
			pdfImp[pos++] = evalPdf(vt, scene, vs, vtPred, EImportance, vtMeasure)
				* sensorSubpath.edge(t - 1)->pdf[EImportance];

			for (int i = t - 1; i>0; --i)
				pdfImp[pos++] = sensorSubpath.vertex(i)->pdf[EImportance]
				* sensorSubpath.edge(i - 1)->pdf[EImportance];
		}

		/* Collect radiance transfer area/volume densities from vertices */
		pos = 0;
		if (s > 0) {
			for (int i = 0; i<s - 1; ++i)
				pdfRad[pos++] = emitterSubpath.vertex(i + 1)->pdf[ERadiance]
				* emitterSubpath.edge(i)->pdf[ERadiance];

			pdfRad[pos++] = evalPdf(vs, scene, vt, vsPred, ERadiance, vsMeasure)
				* emitterSubpath.edge(s - 1)->pdf[ERadiance];
		}

		pdfRad[pos++] = evalPdf(vt, scene, vtPred, vs, ERadiance, vtMeasure)
			* connectionEdge->pdf[ERadiance];

		for (int i = t; i>0; --i)
			pdfRad[pos++] = sensorSubpath.vertex(i - 1)->pdf[ERadiance]
			* sensorSubpath.edge(i - 1)->pdf[ERadiance];

		pdfRad[pos++] = 1.0;


		/* When the path contains specular surface interactions, it is possible
		to compute the correct MI weights even without going through all the
		trouble of computing the proper generalized geometric terms (described
		in the SIGGRAPH 2012 specular manifolds paper). The reason is that these
		all cancel out. But to make sure that that's actually true, we need to
		convert some of the area densities in the 'pdfRad' and 'pdfImp' arrays
		into the projected solid angle measure */
		for (int i = 1; i <= k - 3; ++i) {
			if (i == s || !(connectable[i] && !connectable[i + 1]))
				continue;

			const PathVertex *cur = i <= s ? emitterSubpath.vertex(i) : sensorSubpath.vertex(k - i);
			const PathVertex *succ = i + 1 <= s ? emitterSubpath.vertex(i + 1) : sensorSubpath.vertex(k - i - 1);
			const PathEdge *edge = i < s ? emitterSubpath.edge(i) : sensorSubpath.edge(k - i - 1);

			pdfImp[i + 1] *= edge->length * edge->length / std::abs(
				(succ->isOnSurface() ? dot(edge->d, succ->getGeometricNormal()) : 1) *
				(cur->isOnSurface() ? dot(edge->d, cur->getGeometricNormal()) : 1));
		}

		for (int i = k - 1; i >= 3; --i) {
			if (i - 1 == s || !(connectable[i] && !connectable[i - 1]))
				continue;

			const PathVertex *cur = i <= s ? emitterSubpath.vertex(i) : sensorSubpath.vertex(k - i);
			const PathVertex *succ = i - 1 <= s ? emitterSubpath.vertex(i - 1) : sensorSubpath.vertex(k - i + 1);
			const PathEdge *edge = i <= s ? emitterSubpath.edge(i - 1) : sensorSubpath.edge(k - i);

			pdfRad[i - 1] *= edge->length * edge->length / std::abs(
				(succ->isOnSurface() ? dot(edge->d, succ->getGeometricNormal()) : 1) *
				(cur->isOnSurface() ? dot(edge->d, cur->getGeometricNormal()) : 1));
		}

		int emitterRefIndirection = 2, sensorRefIndirection = k - 2;

		/* One more array sweep before the actual useful work starts -- phew! :)
		"Collapse" edges/vertices that were caused by BSDF::ENull interactions.
		The BDPT implementation is smart enough to connect straight through those,
		so they shouldn't be treated as Dirac delta events in what follows */
		for (int i = 1; i <= k - 3; ++i) {
			if (!connectable[i] || !isNull[i + 1])
				continue;

			int start = i + 1, end = start;
			while (isNull[end + 1])
				++end;

			if (!connectable[end + 1]) {
				/// The chain contains a non-ENull interaction
				isNull[start] = false;
				continue;
			}

			const PathVertex *before = i <= s ? emitterSubpath.vertex(i) : sensorSubpath.vertex(k - i);
			const PathVertex *after = end + 1 <= s ? emitterSubpath.vertex(end + 1) : sensorSubpath.vertex(k - end - 1);

			Vector d = before->getPosition() - after->getPosition();
			Float lengthSquared = d.lengthSquared();
			d /= std::sqrt(lengthSquared);

			Float geoTerm = std::abs(
				(before->isOnSurface() ? dot(before->getGeometricNormal(), d) : 1) *
				(after->isOnSurface() ? dot(after->getGeometricNormal(), d) : 1)) / lengthSquared;

			pdfRad[start - 1] *= pdfRad[end] * geoTerm;
			pdfRad[end] = 1;
			pdfImp[start] *= pdfImp[end + 1] * geoTerm;
			pdfImp[end + 1] = 1;

			/* When an ENull chain starts right after the emitter / before the sensor,
			we must keep track of the reference vertex for direct sampling strategies. */
			if (start == 2)
				emitterRefIndirection = end + 1;
			else if (end == k - 2)
				sensorRefIndirection = start - 1;

			i = end;
		}

		double initial = 1.0f;

		/* When direct sampling strategies are enabled, we must
		account for them here as well */
		if (sampleDirect) {
			/* Direct connection probability of the emitter */
			const PathVertex *sample = s>0 ? emitterSubpath.vertex(1) : vt;
			const PathVertex *ref = emitterRefIndirection <= s
				? emitterSubpath.vertex(emitterRefIndirection) : sensorSubpath.vertex(k - emitterRefIndirection);
			EMeasure measure = sample->getAbstractEmitter()->getDirectMeasure();

			if (connectable[1] && connectable[emitterRefIndirection])
				ratioEmitterDirect = ref->evalPdfDirect(scene, sample, EImportance,
				measure == ESolidAngle ? EArea : measure) / pdfImp[1];

			/* Direct connection probability of the sensor */
			sample = t>0 ? sensorSubpath.vertex(1) : vs;
			ref = sensorRefIndirection <= s ? emitterSubpath.vertex(sensorRefIndirection)
				: sensorSubpath.vertex(k - sensorRefIndirection);
			measure = sample->getAbstractEmitter()->getDirectMeasure();

			if (connectable[k - 1] && connectable[sensorRefIndirection])
				ratioSensorDirect = ref->evalPdfDirect(scene, sample, ERadiance,
				measure == ESolidAngle ? EArea : measure) / pdfRad[k - 1];

			if (s == 1)
				initial /= ratioEmitterDirect;
			else if (t == 1)
				initial /= ratioSensorDirect;
		}

		double weight = 1, pdf = initial;

		/* With all of the above information, the MI weight can now be computed.
		Since the goal is to evaluate the power heuristic, the absolute area
		product density of each strategy is interestingly not required. Instead,
		an incremental scheme can be used that only finds the densities relative
		to the (s,t) strategy, which can be done using a linear sweep. For
		details, refer to the Veach thesis, p.306. */
		for (int i = s + 1; i<k; ++i) {
			double next = pdf * (double)pdfImp[i] / (double)pdfRad[i],
				value = next;

			if (sampleDirect) {
				if (i == 1)
					value *= ratioEmitterDirect;
				else if (i == sensorRefIndirection)
					value *= ratioSensorDirect;
			}


			int tPrime = k - i - 1;
			if (connectable[i] && (connectable[i + 1] || isNull[i + 1]) && (lightImage || tPrime > 1))
				weight += value*value;

			pdf = next;
		}

		/* As above, but now compute pdf[i] with i<s (this is done by
		evaluating the inverse of the previous expressions). */
		pdf = initial;
		for (int i = s - 1; i >= 0; --i) {
			double next = pdf * (double)pdfRad[i + 1] / (double)pdfImp[i + 1],
				value = next;

			if (sampleDirect) {
				if (i == 1)
					value *= ratioEmitterDirect;
				else if (i == sensorRefIndirection)
					value *= ratioSensorDirect;
			}

			int tPrime = k - i - 1;
			if (connectable[i] && (connectable[i + 1] || isNull[i + 1]) && (lightImage || tPrime > 1))
				weight += value*value;

			pdf = next;
		}

		return (Float)(1.0 / weight);
	}

	ref<WorkProcessor> clone() const {
		return new GuidedBDPTRenderer(m_config);
	}

	MTS_DECLARE_CLASS()
private:
	ref<Scene> m_scene;
	ref<Sensor> m_sensor;
	ref<Sampler> m_sampler;
	ref<ReconstructionFilter> m_rfilter;
	MemoryPool m_pool;
	GuidedBDPTConfiguration m_config;
	HilbertCurve2D<uint8_t> m_hilbertCurve;
	ref<GuidingSamplers> m_guidingSampler;
};


/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

GuidedBDPTProcess::GuidedBDPTProcess(const RenderJob *parent, RenderQueue *queue,
	const GuidedBDPTConfiguration &config) :
	BlockedRenderProcess(parent, queue, config.blockSize), m_config(config) {
	m_refreshTimer = new Timer();
}

ref<WorkProcessor> GuidedBDPTProcess::createWorkProcessor() const {
	return new GuidedBDPTRenderer(m_config);
}

void GuidedBDPTProcess::develop() {
	if (!m_config.lightImage)
		return;
	LockGuard lock(m_resultMutex);
	const ImageBlock *lightImage = m_result->getLightImage();
	m_film->setBitmap(m_result->getImageBlock()->getBitmap());
	m_film->addBitmap(lightImage->getBitmap(), 1.0f / m_config.sampleCount);
	m_refreshTimer->reset();
	m_queue->signalRefresh(m_parent);
}

void GuidedBDPTProcess::processResult(const WorkResult *wr, bool cancelled) {
	if (cancelled)
		return;
	const GuidedBDPTWorkResult *result = static_cast<const GuidedBDPTWorkResult *>(wr);
	ImageBlock *block = const_cast<ImageBlock *>(result->getImageBlock());
	LockGuard lock(m_resultMutex);
	m_progress->update(++m_resultCount);
	if (m_config.lightImage) {
		const ImageBlock *lightImage = m_result->getLightImage();
		m_result->put(result);
		if (m_parent->isInteractive()) {
			/* Modify the finished image block so that it includes the light image contributions,
			   which creates a more intuitive preview of the rendering process. This is
			   not 100% correct but doesn't matter, as the shown image will be properly re-developed
			   every 2 seconds and once more when the rendering process finishes */

			Float invSampleCount = 1.0f / m_config.sampleCount;
			const Bitmap *sourceBitmap = lightImage->getBitmap();
			Bitmap *destBitmap = block->getBitmap();
			int borderSize = block->getBorderSize();
			Point2i offset = block->getOffset();
			Vector2i size = block->getSize();

			for (int y=0; y<size.y; ++y) {
				const Float *source = sourceBitmap->getFloatData()
					+ (offset.x + (y+offset.y) * sourceBitmap->getWidth()) * SPECTRUM_SAMPLES;
				Float *dest = destBitmap->getFloatData()
					+ (borderSize + (y + borderSize) * destBitmap->getWidth()) * (SPECTRUM_SAMPLES + 2);

				for (int x=0; x<size.x; ++x) {
					Float weight = dest[SPECTRUM_SAMPLES + 1] * invSampleCount;
					for (int k=0; k<SPECTRUM_SAMPLES; ++k)
						*dest++ += *source++ * weight;
					dest += 2;
				}
			}
		}
	}

	m_film->put(block);

	/* Re-develop the entire image every two seconds if partial results are
	   visible (e.g. in a graphical user interface). This only applies when
	   there is a light image. */
	bool developFilm = m_config.lightImage &&
		(m_parent->isInteractive() && m_refreshTimer->getMilliseconds() > 2000);

	m_queue->signalWorkEnd(m_parent, result->getImageBlock(), false);

	if (developFilm)
		develop();
}

void GuidedBDPTProcess::bindResource(const std::string &name, int id) {
	BlockedRenderProcess::bindResource(name, id);
	if (name == "sensor" && m_config.lightImage) {
		/* If needed, allocate memory for the light image */
		m_result = new GuidedBDPTWorkResult(m_config, NULL, m_film->getCropSize());
		m_result->clear();
	}
}

MTS_IMPLEMENT_CLASS_S(GuidedBDPTRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(GuidedBDPTProcess, false, BlockedRenderProcess)
MTS_NAMESPACE_END
