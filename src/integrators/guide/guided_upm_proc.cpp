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

#include <mitsuba/core/sfcurve.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/bidir/util.h>
#include <mitsuba/bidir/path.h>
#include "guided_upm_proc.h"

#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/fstream.h>

#define EXCLUDE_DIRECT_LIGHTING

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Worker implementation                        */
/* ==================================================================== */

const bool EnableCovAwareMis = true;
const bool MaxClampedConnectionPdf = false;

class GuidedUPMRenderer : public WorkProcessor {
public:
	GuidedUPMRenderer(const GuidedUPMConfiguration &conf)
		: m_config(conf) {
	}

	GuidedUPMRenderer(Stream *stream, InstanceManager *manager)
		: WorkProcessor(stream, manager) {
		m_config = GuidedUPMConfiguration(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		m_config.serialize(stream);
	}

	ref<WorkUnit> createWorkUnit() const {
		return new SeedWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new UPMWorkResult(m_film->getCropSize().x, m_film->getCropSize().y, m_config.maxDepth, m_film->getReconstructionFilter(), true);
	}

	void prepare() {
		Scene *scene = static_cast<Scene *>(getResource("scene"));
		m_sampler = static_cast<Sampler *>(getResource("sampler"));
		m_sensor = static_cast<Sensor *>(getResource("sensor"));
		m_scene = new Scene(scene);
		m_film = m_sensor->getFilm();
		m_scene->setSensor(m_sensor);
		m_scene->setSampler(m_sampler);
		m_scene->removeSensor(scene->getSensor());
		m_scene->addSensor(m_sensor);
		m_scene->setSensor(m_sensor);
		m_scene->wakeup(NULL, m_resources);
		m_scene->initializeBidirectional();

		m_pathSampler = new PathSampler(PathSampler::EBidirectional, m_scene,
			m_sampler, m_sampler, m_sampler, m_config.maxDepth,
			m_config.rrDepth, false /*m_config.separateDirect*/, true /*m_config.directSampling*/,
			true, m_sampler);

		m_guidingSampler = static_cast<GuidingSamplers *>(getResource("guidingSampler"))->clone();
		SLog(EInfo, "clone %d", (void*)m_guidingSampler);
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {
		UPMWorkResult *wr = static_cast<UPMWorkResult *>(workResult);
		wr->clear();
		ImageBlock *midres = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize(), m_film->getReconstructionFilter());
		midres->clear();
		const SeedWorkUnit *wu = static_cast<const SeedWorkUnit *>(workUnit);
		const int workID = wu->getID();
		const int numWork = wu->getTotalWorkNum();
		SplatList *splats = new SplatList();
		splats->clear();

		HilbertCurve2D<int> hilbertCurve;
		TVector2<int> filmSize(m_film->getCropSize());
		hilbertCurve.initialize(filmSize);
		uint64_t iteration = workID;

		int splatcnt = 0;

		size_t actualSampleCount;
		float radius = m_config.initialRadius;
		ref<Timer> timer = new Timer();
		for (actualSampleCount = 0; actualSampleCount < m_config.sampleCount || (wu->getTimeout() > 0 && timer->getSeconds() < wu->getTimeout()); actualSampleCount++) {
			if (m_config.initialRadius > 0.0f){
				Float reduceFactor = 1.0 / std::pow((Float)(iteration + 1), (Float)(0.5 * (1 - m_config.radiusAlpha/*radiusAlpha*/)));
				radius = std::max(reduceFactor * m_config.initialRadius, (Float)1e-7);
				iteration += numWork;
			}

			gatherLightPathsUPM(m_config.useVC, m_config.useVM, radius, hilbertCurve.getPointCount(), wr, m_config.rejectionProb);

			for (size_t i = 0; i < hilbertCurve.getPointCount(); ++i) {
				if (stop) break;

				Point2i offset = Point2i(hilbertCurve[i]);
				m_sampler->generate(offset);
				sampleSplatsUPM(wr, radius, offset, i, *splats, 
					m_config.useVC, m_config.useVM, m_config.rejectionProb, m_config.clampThreshold);

				for (size_t k = 0; k < splats->size(); ++k) {
 					Spectrum value = splats->getValue(k);
					wr->putSample(splats->getPosition(k), &value[0]);
 				}
			}			
		}

		Log(EInfo, "Run %d iterations", actualSampleCount);
		wr->accumSampleCount(actualSampleCount);
		
		delete splats;
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

			if (mode == ERadiance)
				m_guidingSampler->getWeightWindow().pathTracing();
			else
				m_guidingSampler->getWeightWindow().lightTracing();
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
			if (mode != ERadiance)
				m_guidingSampler->getWeightWindow().pathTracing();
			else
				m_guidingSampler->getWeightWindow().lightTracing();
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

	int randomWalk(Path &curPath, const Scene *scene, Sampler *sampler,
		int nSteps, int rrStart, ETransportMode mode,
		MemoryPool &pool) {
		/* Determine the relevant edge and vertex to start the random walk */
		PathVertex *curVertex = curPath.vertex(curPath.vertexCount() - 1),
			*predVertex = curPath.vertexCount() < 2 ? NULL :
			curPath.vertex(curPath.vertexCount() - 2);
		PathEdge *predEdge = (curPath.edgeCount() == 0) ? NULL :
			curPath.edge(curPath.edgeCount() - 1);
		Spectrum throughput(1.0f);

		for (int i = 0; i < nSteps || nSteps == -1; ++i) {
			PathVertex *succVertex = pool.allocVertex();
			PathEdge *succEdge = pool.allocEdge();

			if (!sampleNext(curVertex, scene, sampler, predVertex, predEdge, succEdge,
				succVertex, mode, rrStart != -1 && i >= rrStart, &throughput)) {
				pool.release(succVertex);
				pool.release(succEdge);
				return i;
			}

			curPath.append(succEdge, succVertex);

			predVertex = curVertex;
			curVertex = succVertex;
			predEdge = succEdge;
		}

		return nSteps;
	}
	int randomWalkFromPixel(Path &curPath, const Scene *scene, Sampler *sampler,
		int nSteps, const Point2i &pixelPosition, int rrStart, MemoryPool &pool) {

		PathVertex *v1 = pool.allocVertex(), *v2 = pool.allocVertex();
		PathEdge *e0 = pool.allocEdge(), *e1 = pool.allocEdge();

		/* Use a special sampling routine for the first two sensor vertices so that
		the resulting subpath passes through the specified pixel position */
		int t = curPath.vertex(0)->sampleSensor(scene,
			sampler, pixelPosition, e0, v1, e1, v2);

		if (t < 1) {
			pool.release(e0);
			pool.release(v1);
			return 0;
		}

		curPath.append(e0, v1);

		if (t < 2) {
			pool.release(e1);
			pool.release(v2);
			return 1;
		}

		curPath.append(e1, v2);

		PathVertex *predVertex = v1, *curVertex = v2;
		PathEdge *predEdge = e1;
		Spectrum throughput(1.0f);

		for (; t<nSteps || nSteps == -1; ++t) {
			PathVertex *succVertex = pool.allocVertex();
			PathEdge *succEdge = pool.allocEdge();

			if (!sampleNext(curVertex, scene, sampler, predVertex, predEdge, succEdge,
				succVertex, ERadiance, rrStart != -1 && t >= rrStart, &throughput)) {
				pool.release(succVertex);
				pool.release(succEdge);
				return t;
			}

			curPath.append(succEdge, succVertex);

			predVertex = curVertex;
			curVertex = succVertex;
			predEdge = succEdge;
		}

		return nSteps;
	}

	void gatherLightPathsUPM(const bool useVC, const bool useVM,
		const float gatherRadius, const int nsample, UPMWorkResult *wr, Float rejectionProb){
		const Sensor *sensor = m_scene->getSensor();
		m_pathSampler->m_lightPathTree.clear();
		m_pathSampler->m_lightVertices.clear();
		m_pathSampler->m_lightVerticesExt.clear();
		m_pathSampler->m_lightPathEnds.clear();		

		PathVertex vtPred, vt;
		PathEdge vtEdge, connectionEdge;

		Float time = sensor->getShutterOpen();
		size_t m_lightPathNum = nsample;
		m_pathSampler->m_lightPathNum = m_lightPathNum;
		Float etaVCM = (M_PI * gatherRadius * gatherRadius) * m_lightPathNum * (1.f - rejectionProb);
		Float invLightPathNum = 1.f / m_lightPathNum;
		Float misVmWeightFactor = useVM ? MisHeuristic(etaVCM) : 0.f;
		Float misVcWeightFactor = useVC ? MisHeuristic(1.f / etaVCM) : 0.f;
		for (size_t k = 0; k < m_lightPathNum; k++){
			// emitter states
			MisState emitterState, sensorState;
			Spectrum importanceWeight = Spectrum(1.0f);

			/* Initialize the path endpoints */
			m_pathSampler->m_emitterSubpath.initialize(m_scene, time, EImportance, m_pathSampler->m_pool);

			/* Perform random walks from the emitter side */
			randomWalk(m_pathSampler->m_emitterSubpath, m_scene, m_pathSampler->m_lightPathSampler, m_pathSampler->m_emitterDepth,
				m_pathSampler->m_rrDepth, EImportance, m_pathSampler->m_pool);

			PathVertex* vs = m_pathSampler->m_emitterSubpath.vertex(0);
			LightVertex lvertex = LightVertex(vs, NULL, emitterState, Spectrum(1.f));
			m_pathSampler->m_lightVertices.push_back(lvertex);
			LightVertexExt lvertexExt = LightVertexExt(vs, NULL, 0);
			m_pathSampler->m_lightVerticesExt.push_back(lvertexExt);
			for (int s = 1; s < (int)m_pathSampler->m_emitterSubpath.vertexCount(); ++s) {
				PathVertex
					*vsPred3 = m_pathSampler->m_emitterSubpath.vertexOrNull(s - 3),
					*vsPred2 = m_pathSampler->m_emitterSubpath.vertexOrNull(s - 2),
					*vsPred = m_pathSampler->m_emitterSubpath.vertex(s - 1),
					*vs = m_pathSampler->m_emitterSubpath.vertex(s);
				PathEdge
					*esPred = m_pathSampler->m_emitterSubpath.edgeOrNull(s - 2),
					*es = m_pathSampler->m_emitterSubpath.edgeOrNull(s - 1);

				// Compute the throughput of emitter subpath till this vertex
				importanceWeight *= vsPred->weight[EImportance] * vsPred->rrWeight * es->weight[EImportance];

				// update mis helper and path type
				updateMisHelper(s - 1, m_pathSampler->m_emitterSubpath, emitterState, m_scene, gatherRadius, m_lightPathNum, useVC, useVM, EImportance);

				// store light paths												
				//if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon /* don't save backfaced photons */){
				{
					LightVertex lvertex = LightVertex(vs, vsPred, emitterState, importanceWeight);
					m_pathSampler->m_lightVertices.push_back(lvertex);
					LightVertexExt lvertexExt = LightVertexExt(vs, vsPred, s);
					m_pathSampler->m_lightVerticesExt.push_back(lvertexExt);
					if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon){
						LightPathNode lnode(vs->getPosition(), m_pathSampler->m_lightVertices.size() - 1, s);
						m_pathSampler->m_lightPathTree.push_back(lnode);
					}
				}

				// connect to camera
				if (vs->measure != EDiscrete && wr != NULL && useVC){
					Point2 samplePos(0.0f);
					Spectrum value = importanceWeight * vs->sampleDirect(m_scene, m_pathSampler->m_directSampler,
						&vtPred, &vtEdge, &vt, ERadiance);
					if (value.isZero() || !vt.getSamplePosition(vs, samplePos)) continue;

					value *= vs->eval(m_scene, vsPred, &vt, EImportance);
					vs->measure = EArea;
					int interactions = m_pathSampler->m_maxDepth - s;
					if (value.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, &vt, &vtEdge, interactions))
						continue;
					value *= connectionEdge.evalCached(vs, &vt, PathEdge::ETransmittance | PathEdge::ECosineImp);

					Float miWeight = miWeightVC(m_scene, s, 1, emitterState, sensorState,
						vsPred3, vsPred2, vsPred, vs,
						&vt, &vtPred, NULL, NULL,
						gatherRadius, m_lightPathNum, useVC, useVM);

#if GUPM_DEBUG == 1
					wr->putDebugSample(s, 1, samplePos, value * miWeight);
					wr->putDebugSampleM(s, 1, samplePos, value);
#endif

					value *= miWeight;
					if (value.isZero()) continue;

// 					if ((samplePos.x >= 242 && samplePos.x < 243 && samplePos.y >= 76 && samplePos.y < 77)
// 						&& (s == 4) && value[1] > 0.5f){
// 						float fucka = 1.f;
// 						Float miWeight2 = miWeightVC(m_scene, s, 1, emitterState, sensorState,
// 							vsPred3, vsPred2, vsPred, vs,
// 							&vt, &vtPred, NULL, NULL,
// 							gatherRadius, m_lightPathNum, useVC, useVM);
// 					}

					vt.getSamplePosition(vs, samplePos);
					wr->putSample(samplePos, &value[0]);
				}
			}
			m_pathSampler->m_lightPathEnds.push_back(m_pathSampler->m_lightVertices.size());
		}

		// build kdtree
		m_pathSampler->m_lightPathTree.build(true);

		/* Release any used edges and vertices back to the memory pool */
		m_pathSampler->m_emitterSubpath.release(m_pathSampler->m_pool);
	}

	void sampleSplatsUPM(UPMWorkResult *wr,
		const float gatherRadius, const Point2i &offset,
		const size_t cameraPathIndex, SplatList &list, bool useVC, bool useVM,
		Float rejectionProb, size_t clampThreshold) {
		list.clear();

		const Sensor *sensor = m_scene->getSensor();
		size_t nLightPaths = m_pathSampler->m_lightPathNum;
		Float invLightPaths = 1.f / (float)nLightPaths;

		bool repeated = false;

		PathVertex tempSample, tempEndpoint;
		PathEdge tempEdge;
		switch (m_pathSampler->m_technique) {
		case PathSampler::EBidirectional: {
			/* Uniformly sample a scene time */
			Float time = sensor->getShutterOpen();
			if (sensor->needsTimeSample())
				time = sensor->sampleTime(m_pathSampler->m_sensorSampler->next1D());

			/* Initialize the path endpoints */
			m_pathSampler->m_sensorSubpath.initialize(m_scene, time, ERadiance, m_pathSampler->m_pool);

			if (offset == Point2i(-1))
				randomWalk(m_pathSampler->m_sensorSubpath, m_scene, m_pathSampler->m_sensorSampler,
				m_pathSampler->m_sensorDepth, m_pathSampler->m_rrDepth, ERadiance, m_pathSampler->m_pool);
			else
				randomWalkFromPixel(m_pathSampler->m_sensorSubpath, m_scene, m_pathSampler->m_sensorSampler,
				m_pathSampler->m_sensorDepth, offset, m_pathSampler->m_rrDepth, m_pathSampler->m_pool);

			/* Compute the combined weights along the two subpaths */
			Spectrum *radianceWeights = (Spectrum *)alloca(m_pathSampler->m_sensorSubpath.vertexCount()  * sizeof(Spectrum));

			radianceWeights[0] = Spectrum(1.0f);

			for (size_t i = 1; i<m_pathSampler->m_sensorSubpath.vertexCount(); ++i)
				radianceWeights[i] = radianceWeights[i - 1] *
				m_pathSampler->m_sensorSubpath.vertex(i - 1)->weight[ERadiance] *
				m_pathSampler->m_sensorSubpath.vertex(i - 1)->rrWeight *
				m_pathSampler->m_sensorSubpath.edge(i - 1)->weight[ERadiance];

			bool watchThread = false;
			Point2 initialSamplePos(0.0f);
			if (m_pathSampler->m_sensorSubpath.vertexCount() > 2) {
				Point2 samplePos(0.0f);
				m_pathSampler->m_sensorSubpath.vertex(1)->getSamplePosition(m_pathSampler->m_sensorSubpath.vertex(2), samplePos);
				list.append(samplePos, Spectrum(0.0f));
				initialSamplePos = samplePos;
			}

			// initialize of MIS helper
			// VCM tech report Eq. 31 - 33		
			MisState *sensorStates = (MisState *)alloca(m_pathSampler->m_sensorSubpath.vertexCount() * sizeof(MisState));
			initializeMisHelper(m_pathSampler->m_sensorSubpath, sensorStates, m_scene,
				gatherRadius, nLightPaths, useVC, useVM, ERadiance);

			// massive vertex merging
			if (useVM){
				PathVertex *vs_ = m_pathSampler->m_pool.allocVertex();
				PathVertex *vsPred_ = m_pathSampler->m_pool.allocVertex();
				PathVertex *vsPred2_ = m_pathSampler->m_pool.allocVertex();
				PathVertex *vsPred3_ = m_pathSampler->m_pool.allocVertex();
				PathVertex *vs = NULL, *vsPred = NULL, *vsPred2 = NULL, *vsPred3 = NULL;
				PathEdge connectionEdge;
				PathVertex *succVertex = m_pathSampler->m_pool.allocVertex();
				PathEdge *succEdge = m_pathSampler->m_pool.allocEdge();
				Point2 samplePos(0.0f);
				std::vector<uint32_t> searchResults;
				//std::vector<Point> searchPos;
				std::vector<uint32_t> acceptCnt;
				std::vector<size_t> shootCnt;

				int minT = 2; int minS = 2;

				int maxT = (int)m_pathSampler->m_sensorSubpath.vertexCount() - 1;
				if (m_pathSampler->m_maxDepth != -1)
					maxT = std::min(maxT, m_pathSampler->m_maxDepth);
				for (int t = minT; t <= maxT; ++t) {
					PathVertex
						*vtPred3 = m_pathSampler->m_sensorSubpath.vertexOrNull(t - 3),
						*vtPred2 = m_pathSampler->m_sensorSubpath.vertex(t - 2),
						*vtPred = m_pathSampler->m_sensorSubpath.vertex(t - 1),
						*vt = m_pathSampler->m_sensorSubpath.vertex(t);
					PathEdge
						*predEdge = m_pathSampler->m_sensorSubpath.edge(t - 2);

					if (!vt->isDegenerate()){
						BDAssert(vt->type == PathVertex::ESurfaceInteraction);

						searchResults.clear();
						m_pathSampler->m_lightPathTree.search(vt->getPosition(), gatherRadius, searchResults);						
						if (searchResults.size() == 0){
							continue;
						}

						//searchPos.resize(searchResults.size());
						acceptCnt.resize(searchResults.size());
						shootCnt.resize(searchResults.size());
						for (int i = 0; i < searchResults.size(); i++){
							acceptCnt[i] = 0;
							shootCnt[i] = 0;
						}

						// evaluate sampling domain pdf normalization
						int shareShootThreshold = 32;
						Float invBrdfIntegral = 1.f;
						bool shareShoot = false;			

						MisState sensorState = sensorStates[t - 1];
						MisState sensorStatePred = sensorStates[t - 2];
						for (int k = 0; k < searchResults.size(); k++){
							LightPathNode node = m_pathSampler->m_lightPathTree[searchResults[k]];
							int s = node.data.depth;
							if (m_pathSampler->m_maxDepth != -1 && s + t > m_pathSampler->m_maxDepth + 2 || s < minS) continue;

#ifdef EXCLUDE_DIRECT_LIGHTING
							if (s == 2 && t == 2) continue;
#endif

							if (m_pathSampler->m_sensorSampler->next1D() < rejectionProb) continue;

							size_t vertexIndex = node.data.vertexIndex;
							LightVertex vi = m_pathSampler->m_lightVertices[vertexIndex];
							LightVertex viPred = m_pathSampler->m_lightVertices[vertexIndex - 1];
							MisState emitterState = vi.emitterState;
							MisState emitterStatePred = viPred.emitterState;

							vs = vs_; vsPred = vsPred_; vsPred2 = vsPred2_; vsPred3 = vsPred3_;
							m_pathSampler->m_lightVerticesExt[vertexIndex].expand(vs);
							m_pathSampler->m_lightVerticesExt[vertexIndex - 1].expand(vsPred);
							if (s > 2){
								m_pathSampler->m_lightVerticesExt[vertexIndex - 2].expand(vsPred2);
								m_pathSampler->m_lightVerticesExt[vertexIndex - 3].expand(vsPred3);
							}
							else if (s == 2){
								m_pathSampler->m_lightVerticesExt[vertexIndex - 2].expand(vsPred2);
								vsPred3 = NULL;
							}
							else{
								vsPred2 = NULL;
								vsPred3 = NULL;
							}

							// decide the direction to do connection
							bool cameraDirConnection = (connectionDirection(vsPred, vtPred) == ERadiance);

							samplePos = initialSamplePos;
							if (vtPred->isSensorSample()){
								if (!vtPred->getSamplePosition(cameraDirConnection ? vs : vt, samplePos))
									continue;
							}

							// evaluate contribution
							Spectrum contrib;
							if (cameraDirConnection){
								contrib = radianceWeights[t - 1] * vi.importanceWeight * invLightPaths;
								contrib *= vs->eval(m_scene, vsPred, vtPred, EImportance) *	vtPred->eval(m_scene, vtPred2, vs, ERadiance);
								int interactions = m_pathSampler->m_maxDepth - s - t + 1;
								if (contrib.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, vtPred, NULL, interactions))
									continue;
								contrib *= connectionEdge.evalCached(vs, vtPred, PathEdge::EGeneralizedGeometricTerm);
							}
							else{
								contrib = radianceWeights[t] * viPred.importanceWeight * invLightPaths;
								contrib *= vt->eval(m_scene, vtPred, vsPred, ERadiance) *	vsPred->eval(m_scene, vsPred2, vt, EImportance);
								int interactions = m_pathSampler->m_maxDepth - s - t + 1;
								if (contrib.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vt, vsPred, NULL, interactions))
									continue;
								contrib *= connectionEdge.evalCached(vt, vsPred, PathEdge::EGeneralizedGeometricTerm);
							}

							contrib /= (1.f - rejectionProb);

							Float invp = 0.f;
							if (shareShoot){
								invp = (acceptCnt[k] > 0) ? (Float)(shootCnt[k]) / (Float)(acceptCnt[k]) * invBrdfIntegral : 0;
								contrib *= invp;
							}
							else{
								std::vector<Vector2> componentCDFs;
								std::vector<Vector4> componentBounds;
								Float brdfIntegral;
								
								bool gsampler_valid = (cameraDirConnection && vtPred->isSurfaceInteraction() || !cameraDirConnection && vsPred->isSurfaceInteraction());
								Intersection gits = (cameraDirConnection) ? vtPred->getIntersection() : vsPred->getIntersection();
								GuidedBRDF gsampler(gits, (cameraDirConnection) ? m_guidingSampler->getRadianceSampler() : m_guidingSampler->getImportanceSampler(),
									m_guidingSampler->getConfig().m_mitsuba.bsdfSamplingProbability, gsampler_valid);

								if (cameraDirConnection)
									brdfIntegral = gatherAreaPdf(vtPred, vs->getPosition(), gatherRadius, vtPred2, componentCDFs, componentBounds, gsampler);
								else
									brdfIntegral = gatherAreaPdf(vsPred, vt->getPosition(), gatherRadius, vsPred2, componentCDFs, componentBounds, gsampler);

								if (brdfIntegral == 0.f) continue;
								invBrdfIntegral = 1.f / brdfIntegral;
								size_t totalShoot = 0, acceptedShoot = 0, targetShoot = 1;
								Float distSquared = gatherRadius * gatherRadius;
								while (totalShoot < clampThreshold){
									totalShoot++;

									// restricted sampling evaluation shoots
									Float pointDistSquared;
									if (cameraDirConnection){
										if (!sampleShoot(vtPred, m_scene, m_pathSampler->m_sensorSampler, vtPred2, predEdge, succEdge, succVertex, ERadiance, vs->getPosition(), gatherRadius, componentCDFs, componentBounds, gsampler))
											continue;
										pointDistSquared = (succVertex->getPosition() - vs->getPosition()).lengthSquared();
									}
									else{
										if (!sampleShoot(vsPred, m_scene, m_pathSampler->m_emitterSampler, vsPred2, predEdge, succEdge, succVertex, EImportance, vt->getPosition(), gatherRadius, componentCDFs, componentBounds, gsampler))
											continue;
										pointDistSquared = (succVertex->getPosition() - vt->getPosition()).lengthSquared();
									}

									if (pointDistSquared < distSquared){
										acceptedShoot++;
										if (acceptedShoot == targetShoot){
											break;
										}
									}
								}
								invp = (acceptedShoot > 0) ? (Float)(totalShoot) / (Float)(acceptedShoot)* invBrdfIntegral : 0;
								contrib *= invp;
							}

							// accumulate to image
							if (contrib.isZero()) continue;

							// MIS weighting						
							Float miWeight = miWeightVM(m_scene, s, t,
								emitterState, sensorState, emitterStatePred, sensorStatePred,
								vsPred3, vsPred2, vsPred, vs,
								vt, vtPred, vtPred2, vtPred3,
								cameraDirConnection, gatherRadius, m_pathSampler->m_lightPathNum, useVC, useVM);

#if GUPM_DEBUG == 1
							if (cameraDirConnection){
								wr->putDebugSample(s, t - 1, samplePos, contrib * miWeight);
								wr->putDebugSampleM(s, t - 1, samplePos, contrib);
							}
							else{
								wr->putDebugSample(s - 1, t, samplePos, contrib * miWeight);
								wr->putDebugSampleM(s - 1, t, samplePos, contrib);
							}
#endif

							contrib *= miWeight;

#ifdef GUPM_DEBUG_HARD
							if (contrib[0] < 0.f || _isnan(contrib[0]) || contrib[0] > 100000000.f){
								Log(EWarn, "Invalid sample value[UPM]: %f %f %f, invp = %f, miWeight = %f",
									contrib[0], contrib[1], contrib[2], invp, miWeight);
								continue;
							}
#endif

							if (t == 2) {
								list.append(samplePos, contrib);
							}
							else {
								list.accum(0, contrib);
							}
						}
					}
				}
				m_pathSampler->m_pool.release(vs_);
				m_pathSampler->m_pool.release(vsPred_);
				m_pathSampler->m_pool.release(vsPred2_);
				m_pathSampler->m_pool.release(vsPred3_);
				m_pathSampler->m_pool.release(succVertex);
				m_pathSampler->m_pool.release(succEdge);
			}

			repeated = false;

			// vertex connection	 
			if (useVC){
				PathVertex *vs_ = m_pathSampler->m_pool.allocVertex();
				PathVertex *vsPred_ = m_pathSampler->m_pool.allocVertex();
				PathVertex *vsPred2_ = m_pathSampler->m_pool.allocVertex();
				PathVertex *vsPred3_ = m_pathSampler->m_pool.allocVertex();
				size_t lightPathBegin = (cameraPathIndex == 0) ? 0 : m_pathSampler->m_lightPathEnds[cameraPathIndex - 1];
				size_t lightPathEnd = m_pathSampler->m_lightPathEnds[cameraPathIndex];
				Point2 samplePos(0.0f);
				PathEdge connectionEdge;
				for (size_t i = lightPathBegin; i < lightPathEnd + 2; i++){
					int s = lightPathEnd + 1 - i;
					PathVertex *vsPred3 = vsPred3_, *vsPred2 = vsPred2_, *vsPred = vsPred_, *vs = vs_;
					Spectrum importanceWeight = Spectrum(1.0f);
					MisState emitterState;

					memset(vsPred, 0, sizeof(PathVertex));
					memset(vs, 0, sizeof(PathVertex));
					if (i < lightPathEnd){
						LightVertex lvertex = m_pathSampler->m_lightVertices[i];
						s = m_pathSampler->m_lightVerticesExt[i].depth;
						if (s <= 1) continue;

						importanceWeight = lvertex.importanceWeight;
						emitterState = lvertex.emitterState;

						m_pathSampler->m_lightVerticesExt[i].expand(vs);
						m_pathSampler->m_lightVerticesExt[i - 1].expand(vsPred);
						if (s > 2){
							m_pathSampler->m_lightVerticesExt[i - 2].expand(vsPred2);
							m_pathSampler->m_lightVerticesExt[i - 3].expand(vsPred3);
						}
						else if (s == 2){
							m_pathSampler->m_lightVerticesExt[i - 2].expand(vsPred2);
							vsPred3 = NULL;
						}
						else{
							vsPred2 = NULL;
							vsPred3 = NULL;
						}
					}

					int minT = 2, maxT = (int)m_pathSampler->m_sensorSubpath.vertexCount() - 1;
					if (m_pathSampler->m_maxDepth != -1)
						maxT = std::min(maxT, m_pathSampler->m_maxDepth + 1 - s);
					for (int t = minT; t <= maxT; ++t) {
						PathVertex
							*vtPred3 = m_pathSampler->m_sensorSubpath.vertexOrNull(t - 3),
							*vtPred2 = m_pathSampler->m_sensorSubpath.vertex(t - 2),
							*vtPred = m_pathSampler->m_sensorSubpath.vertex(t - 1),
							*vt = m_pathSampler->m_sensorSubpath.vertex(t);
						PathEdge
							*vtEdge = m_pathSampler->m_sensorSubpath.edgeOrNull(t - 1);

						/* Will receive the path weight of the (s, t)-connection */
						Spectrum value = Spectrum(0.0f);

						RestoreMeasureHelper rmh1(vt);

						/* Will be set to true if direct sampling was used */
						bool sampleDirect = false;

						/* Number of edges of the combined subpaths */
						int depth = s + t - 1;

						/* Allowed remaining number of ENull vertices that can
						be bridged via pathConnect (negative=arbitrarily many) */
						int remaining = m_pathSampler->m_maxDepth - depth;

						/* Account for the terms of the measurement contribution
						function that are coupled to the connection endpoints */
						if (s == 0){
							/* If possible, convert 'vt' into an emitter sample */
							if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
								continue;
							vs = vs_;
							vs->makeEndpoint(m_scene, time, EImportance);
							value = radianceWeights[t] *
								vs->eval(m_scene, NULL, vt, EImportance) *
								vt->eval(m_scene, vtPred, vs, ERadiance);
						}
						else if (s == 1) {
							if (vt->isDegenerate())
								continue;
							/* Generate a position on an emitter using direct sampling */
							value = radianceWeights[t] * vt->sampleDirect(m_scene, m_pathSampler->m_directSampler,
								&tempEndpoint, &tempEdge, &tempSample, EImportance);
							if (value.isZero())
								continue;
							vs = &tempSample; vsPred = &tempEndpoint;
							value *= vt->eval(m_scene, vtPred, vs, ERadiance);
							vt->measure = EArea;
							sampleDirect = true;
						}
						else {
							/* Can't connect degenerate endpoints */
							if ((vs->isDegenerate()) || vt->isDegenerate())
								continue;

							value = importanceWeight * radianceWeights[t] *
								vs->eval(m_scene, vsPred, vt, EImportance) *
								vt->eval(m_scene, vtPred, vs, ERadiance);

							/* Temporarily force vertex measure to EArea. Needed to
							handle BSDFs with diffuse + specular components */
							vs->measure = vt->measure = EArea;
						}

						/* Attempt to connect the two endpoints, which could result in
						the creation of additional vertices (index-matched boundaries etc.) */
						int interactions = remaining;
						if (value.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, vt, NULL, interactions))
							continue;

						/* Determine the pixel sample position when necessary */
						samplePos = initialSamplePos;

						/* Account for the terms of the measurement contribution
						function that are coupled to the connection edge */
						if (!sampleDirect)
							value *= connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
						else
							value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
							(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

						MisState sensorState = sensorStates[t - 1];
						Float miWeight = miWeightVC(m_scene, s, t, emitterState, sensorState,
							vsPred3, vsPred2, vsPred, vs,
							vt, vtPred, vtPred2, vtPred3,
							gatherRadius, m_pathSampler->m_lightPathNum, useVC, useVM);

#if GUPM_DEBUG == 1
						wr->putDebugSample(s, t, samplePos, value * miWeight);
						wr->putDebugSampleM(s, t, samplePos, value);
#endif			

						value *= miWeight;
						if (value.isZero()) continue;



#ifdef GUPM_DEBUG_HARD
						if (value[0] < 0.f || _isnan(value[0]) || value[0] > 100000000.f){
							Log(EWarn, "Invalid sample value[VC]: %f %f %f", value[0], value[1], value[2]);
							continue;
						}
#endif

						if (t < 2) {
							list.append(samplePos, value);
						}
						else {
							list.accum(0, value);
						}
					}
				}
				m_pathSampler->m_pool.release(vs_);
				m_pathSampler->m_pool.release(vsPred_);
				m_pathSampler->m_pool.release(vsPred2_);
				m_pathSampler->m_pool.release(vsPred3_);
			}

			/* Release any used edges and vertices back to the memory pool */
			m_pathSampler->m_sensorSubpath.release(m_pathSampler->m_pool);
		}
			break;

		default:
			Log(EError, "PathSampler::sample(): invalid technique!");
		}
	}

	Float gatherAreaPdf(PathVertex* current, Point p, Float radius, PathVertex* pPred,
		std::vector<Vector2> &componentCDFs, std::vector<Vector4> &componentBounds,
		GuidedBRDF gsampler){
		switch (current->type) {
		case PathVertex::ESensorSample: {
			// Assume perspective camera
			PositionSamplingRecord &pRec = current->getPositionSamplingRecord();
			const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
			Vector4 bbox = sensor->evaluateSphereBounds(p, radius);
			Float prob = (bbox.y - bbox.x) * (bbox.w - bbox.z);
			componentCDFs.push_back(Vector2(prob, 0.f));
			componentBounds.push_back(bbox);
			return prob;
		}
		case PathVertex::ESurfaceInteraction: {
// 			const Intersection &its = current->getIntersection();
// 			Vector wo = p - its.p;
// 			Vector wi = normalize(pPred->getPosition() - its.p);
// 			const BSDF *bsdf = its.getBSDF();
// 			return bsdf->gatherAreaPdf(its.toLocal(wi), its.toLocal(wo), radius, componentCDFs, componentBounds);

			if (true){
				Vector2 rootnode = Vector2();
				int numNode = 2;
				componentCDFs.push_back(Vector2(0.f/* toal pdf, later fill in */, *(float*)&numNode));		// level root node			
				componentCDFs.push_back(Vector2(0.f/* bsdf pdf, later fill in */, 0.f/* pointer to bsdf node, later fill in */));			// bsdf node
				componentCDFs.push_back(Vector2(0.f/* GMM pdf, later fill in */, 0.f/* pointer to GMM node, later fill in */));			// GMM node
				Float bsdfSamplingWeight = m_guidingSampler->getConfig().m_mitsuba.bsdfSamplingProbability;

				if (gsampler.m_pureSpecularBSDF || gsampler.m_impDistrib == NULL)
					bsdfSamplingWeight = 1.f;

				int ptrNode = componentCDFs.size();
				componentCDFs[1].y = *(float*)&ptrNode;
				const Intersection &its = current->getIntersection();
				Vector wo = p - its.p;
				Vector wi = normalize(pPred->getPosition() - its.p);
				const BSDF *bsdf = its.getBSDF();
				Float probBsdf = bsdf->gatherAreaPdf(its.toLocal(wi), its.toLocal(wo), radius, componentCDFs, componentBounds);
				componentCDFs[1].x = probBsdf * bsdfSamplingWeight;

				ptrNode = componentCDFs.size();
				componentCDFs[2].y = *(float*)&ptrNode;
				Float probGMM = 1.f;
				componentCDFs[2].x = probGMM * (1.f - bsdfSamplingWeight);

				Float prob = probBsdf * bsdfSamplingWeight + probGMM * (1.f - bsdfSamplingWeight);
				componentCDFs[0].x = prob;
				return prob;
			}
		}
		case PathVertex::EEmitterSample: {
			// assume no sampling from emitter, bounded CDF and bounded sampling left untouched
			BDAssert(false);
			// 		PositionSamplingRecord &pRec = getPositionSamplingRecord();
			// 		const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			// 		return emitter->gatherAreaPdf(pRec, p, radius, componentCDFs, componentBounds);
		}
		default:
			SLog(EError, "PathVertex::sampsamplingProbabilityleNext(): Encountered an "
				"unsupported vertex type (%i)!", current->type);
			return 0.f;
		}
	}

	bool sampleShoot(PathVertex* current, const Scene *scene, Sampler *sampler,
		const PathVertex *pred, const PathEdge *predEdge,
		PathEdge *succEdge, PathVertex *succ,
		ETransportMode mode, Point gatherPosition, Float gatherRadius,
		std::vector<Vector2> componentCDFs, std::vector<Vector4> componentBounds, 
		GuidedBRDF gsampler) {
		Ray ray;

		memset(succEdge, 0, sizeof(PathEdge));
		memset(succ, 0, sizeof(PathVertex));

		size_t totalSmpl = 0, clampThreshold = 10000000;

		switch (current->type) {
		case PathVertex::ESensorSample: {
			BDAssert(mode == ERadiance && pred->type == PathVertex::ESensorSupernode);
			PositionSamplingRecord &pRec = current->getPositionSamplingRecord();
			const Sensor *sensor = static_cast<const Sensor *>(pRec.object);
			DirectionSamplingRecord dRec;

			/* Sample the image plane */
			Point2 smp = sampler->next2D();
			Vector4 bbox = componentBounds[0];
			smp.x = (bbox.y - bbox.x) * smp.x + bbox.x;
			smp.y = (bbox.w - bbox.z) * smp.y + bbox.z;
			Spectrum result = sensor->sampleDirection(dRec, pRec, smp);

			ray.time = pRec.time;
			ray.setOrigin(pRec.p);
			ray.setDirection(dRec.d);
		}
			break;

		case PathVertex::ESurfaceInteraction: {
			const Intersection &its = current->getIntersection();
			const BSDF *bsdf = its.getBSDF();
			Vector wi = normalize(pred->getPosition() - its.p);
			Vector wo = gatherPosition - its.p;

			// sample CDF tree
			Float rndsmp = sampler->next1D();
			Vector2 rootnode = componentCDFs[0];
			Float invTotalPdf = 1.f / rootnode.x;
			int numNode = *(int*)&rootnode.y;
			Float cdfi = 0.f;
			int chosenLobe = -1;
			int ptrNode = -1;
			for (int i = 0; i < numNode; i++){
				Vector2 nodei = componentCDFs[1 + i];
				Float pdfi = nodei.x;
				if (rndsmp <= (cdfi + pdfi) * invTotalPdf){
					// choose this component
					chosenLobe = i;
					ptrNode = *(int*)&nodei.y;
					//rndsmp = (rndsmp - cdfi * invTotalPdf) / (pdfi * invTotalPdf);
					break;
				}
				cdfi += pdfi;
			}
			if (chosenLobe == 0){
				/* Sample the BSDF */
				Vector dir = bsdf->sampleGatherArea(its.toLocal(wi), its.toLocal(wo), gatherRadius, sampler->next2D(),
					ptrNode, componentCDFs, componentBounds);
				if (dir == Vector(0.f)) return false;
				wo = its.toWorld(dir);				
			}
			else{
				Float pdf = 0.f;
				gsampler.sampleGMM(wo, pdf, sampler);
				if (pdf == 0.f) return false;
			}
			ray.time = its.time;
			ray.setOrigin(its.p);
			ray.setDirection(wo);
		}
			break;

		case PathVertex::EEmitterSample: {
			// assume no sampling from emitter, bounded CDF and bounded sampling left untouched
			BDAssert(false);
			// 		PositionSamplingRecord &pRec = getPositionSamplingRecord();
			// 		const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			// 		DirectionSamplingRecord dRec;
			// 
			// 		Vector dir = emitter->sampleGatherArea(dRec, pRec, gatherPosition, gatherRadius, sampler->next2D(), componentProbs, componentBounds);
			// 		if (dir == Vector(0.f)) return false;
			// 
			// 		ray.time = pRec.time;
			// 		ray.setOrigin(pRec.p);
			// 		ray.setDirection(dir);
		}
			break;

		default:
			SLog(EError, "PathVertex::sampleNext(): Encountered an "
				"unsupported vertex type (%i)!", current->type);
			return false;
		}

		if (!succEdge->sampleNext(scene, sampler, current, ray, succ, mode)) {
			/* Sampling a successor edge + vertex failed, hence the vertex
			is not committed to a particular measure yet -- revert. */
			return false;
		}

		return true;
	}

	Float evalPdf(const PathVertex* current, const Scene *scene, const PathVertex *pred,
		const PathVertex *succ, ETransportMode mode, EMeasure measure) {

		//return current->evalPdf(scene, pred, succ, mode, measure);

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
			if (mode == ERadiance)
				m_guidingSampler->getWeightWindow().pathTracing();
			else
				m_guidingSampler->getWeightWindow().lightTracing();
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

	inline Float MisHeuristic(Float pdf) {
		return pdf * pdf;
	}
	ETransportMode connectionDirection(const PathVertex* vsPred, const PathVertex* vtPred){
		// true for camera direction and false for light direction
		bool cameraDirConnection = true;
		Float sBandwidth = 0.f, tBandwidth = 0.f;
		if (vtPred->isSensorSample()){
			tBandwidth = 10000.f;
		}
		else if (vtPred->isSurfaceInteraction()){
			const Intersection &its = vtPred->getIntersection();
			const BSDF *bsdf = its.getBSDF();
			tBandwidth = bsdf->getBandwidth();
		}
		else
			return ETransportModes;

		if (vsPred->isEmitterSample()){
			const PositionSamplingRecord &pRec = vsPred->getPositionSamplingRecord();
			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);

			sBandwidth = emitter->getBandwidth();
		}
		else if (vsPred->isSurfaceInteraction()){
			const Intersection &its = vsPred->getIntersection();
			const BSDF *bsdf = its.getBSDF();
			sBandwidth = bsdf->getBandwidth();
		}
		else
			return ETransportModes;

		if (sBandwidth < tBandwidth)
			cameraDirConnection = false;

		return cameraDirConnection ? ERadiance : EImportance;
	}
	
	Float misEffectiveEta(int i, Float pi, Float pir, const PathVertex* vPred, const PathVertex* vNext,
		Float gatherRadius, size_t numLightPath, ETransportMode mode, int j = 999){

		if (i < 1 || j < 1) return 0.f;
#ifdef EXCLUDE_DIRECT_LIGHTING
		if (i + j <= 2) return 0.f; // exclude (2,2) photon mapping
#endif

#ifdef EXCLUDE_DIRECT_LIGHT_UPM
		if (i == 1 && mode == EImportance)
			return 0.f;
#endif

		Float eta = 0.f;
		ETransportMode connectionMode = (mode == EImportance) ? connectionDirection(vPred, vNext) : connectionDirection(vNext, vPred);
		Float ps = (connectionMode != mode) ? pir : pi;
		if (EnableCovAwareMis){
			//Float covfactor = (pc == 0.f) ? 0.f : std::min(1.f, 1.f / (pc * M_PI * gatherRadius * gatherRadius));
			Float covfactor = 1.f;
			eta = (ps == 0.f) ? 0.f : std::min(1.f, ps * M_PI * gatherRadius * gatherRadius) * ((Float)numLightPath * covfactor) / ps;
		}
		else
			eta = M_PI * gatherRadius * gatherRadius * (Float)numLightPath;
		return eta;
	}
	bool connectable(const PathVertex* v, const PathVertex* vNext = NULL){
		if (v->isEmitterSupernode() || v->isSensorSupernode()){
			BDAssert(vNext != NULL);
			const AbstractEmitter *emitter = vNext->getAbstractEmitter();
			EMeasure emitterDirectMeasure = emitter->getDirectMeasure();
			return emitterDirectMeasure != EDiscrete && emitterDirectMeasure != EInvalidMeasure;
		}
		else if (v->isEmitterSample() || v->isSensorSample()){
			const AbstractEmitter *emitter = v->getAbstractEmitter();
			EMeasure emitterDirectMeasure = emitter->getDirectMeasure();
			return emitterDirectMeasure != EInvalidMeasure;
		}
		return v->isConnectable();
	}
	Float misEVC(int i, MisState statePred, Float piPred, Float pir2, const PathVertex* vPred, const PathVertex* vPred2){

		Float eVC = statePred[EVC];
		if (i == 0)
			return 0.f;
		else{
			Float ratioDirect = 1.f;
			if (i == 2)
				ratioDirect = statePred[EDIR];
			Float invPiPred = (piPred == 0.f) ? 0.f : (1.f / piPred);
			eVC = MisHeuristic(pir2 * invPiPred) * eVC;
			if (connectable(vPred) && connectable(vPred2, vPred)){ // TO TEST: camera vertex, point light, dir light, area light
				eVC += MisHeuristic(invPiPred) * ratioDirect;
			}
		}
		return eVC;
	}
	Float misEPM(int i, MisState statePred, Float piPred, Float pir2, const PathVertex* vPred, const PathVertex* vPred3,
		Float gatherRadius, size_t numLightPath, ETransportMode mode){

		Float ePM = 0.f;
		if (i < 3)
			return 0.f;
#ifdef EXCLUDE_DIRECT_LIGHT_UPM
		else if (i == 3)
			return 0.f;
#endif
		else{
			ePM = statePred[EPM];
			Float invPiPred = (piPred == 0.f) ? 0.f : (1.f / piPred);
			ePM *= MisHeuristic(pir2 * invPiPred);

			Float piPred2 = vPred3->pdf[mode];
			Float etar2 = misEffectiveEta(i - 2, piPred2, pir2, vPred3, vPred, gatherRadius, numLightPath, mode);
			ePM += MisHeuristic(pir2 * invPiPred * etar2);
		}
		return ePM;
	}
	Float misWeightVC_vc(int i, MisState statei,
		Float pi, Float pir, Float pir1,
		const PathVertex* v, const PathVertex* vPred){

		Float invpi = (pi == 0.f) ? 0.f : 1.f / pi;
		Float ratioDirect = (i == 1) ? statei[EDIR] : 1.f;

		Float wVC = statei[EVC] * MisHeuristic(pir * invpi * pir1);
		if (connectable(v) && connectable(vPred, v)){ // TO TEST: camera vertex, point light, dir light, area light
			wVC += MisHeuristic(pir * invpi) * ratioDirect;
		}
		return wVC;
	}
	Float misWeightVC_pm(int i, int j, MisState statei,
		Float pi, Float piPred, Float pir, Float pir1,
		const PathVertex* vNext, const PathVertex* v, const PathVertex* vPred, const PathVertex* vPred2,
		Float gatherRadius, size_t numLightPath, ETransportMode mode){

		Float wPM = 0.f;

		Float etar1 = misEffectiveEta(i - 1, piPred, pir1, vPred2, v, gatherRadius, numLightPath, mode, j + 2);
		Float invpi = (pi == 0.f) ? 0.f : 1.f / pi;
		wPM = MisHeuristic(pir * invpi) * (MisHeuristic(pir1) * statei[EPM] + MisHeuristic(pir1 * etar1));

		Float etar = misEffectiveEta(i, pi, pir, vPred, vNext, gatherRadius, numLightPath, mode, j + 1);
		wPM += MisHeuristic(pir * etar);

		return wPM;
	}
	Float misWeightPM(int i, int j, MisState statei,
		Float pi, Float piPred, Float pir, Float pir1,
		const PathVertex* vNext, const PathVertex* v, const PathVertex* vPred, const PathVertex* vPred2,
		Float gatherRadius, size_t numLightPath, ETransportMode mode,
		bool useVC, bool useVM){

		Float etar = misEffectiveEta(i, pi, pir, vPred, vNext, gatherRadius, numLightPath, mode, j);
		Float etar1 = misEffectiveEta(i - 1, piPred, pir1, vPred2, v, gatherRadius, numLightPath, mode, j + 1);
		Float invpi = (pi == 0.f) ? 0.f : 1.f / pi;
		Float invetar = (etar == 0.f) ? 0.f : 1.f / etar;
		Float ratioDirect = (i == 1) ? statei[EDIR] : 1.f;

		Float wPM = MisHeuristic(invpi * invetar) * (MisHeuristic(pir1 * etar1) + MisHeuristic(pir1) * statei[EPM]);
		Float wVC = MisHeuristic(invpi * invetar) * (ratioDirect + MisHeuristic(pir1) * statei[EVC]);

		return (useVC ? wVC : 0.f) + (useVM ? wPM : 0.f);
	}
	Float misWeightPM_pred(int i, int j, MisState state, MisState statePred,
		Float pi, Float piPred, Float piPred2, Float pir, Float pir1, Float pir2,
		const PathVertex* vNext, const PathVertex* v, const PathVertex* vPred, const PathVertex* vPred2, const PathVertex* vPred3,
		Float gatherRadius, size_t numLightPath, ETransportMode mode,
		bool useVC, bool useVM){

		Float etar = misEffectiveEta(i, pi, pir, vPred, vNext, gatherRadius, numLightPath, mode, j);
		Float etar1 = misEffectiveEta(i - 1, piPred, pir1, vPred2, v, gatherRadius, numLightPath, mode, j + 1);
		Float etar2 = misEffectiveEta(i - 2, piPred2, pir2, vPred3, vPred, gatherRadius, numLightPath, mode, j + 2);
		Float invpi = (pi == 0.f) ? 0.f : 1.f / pi;
		Float invpiPred = (piPred == 0.f) ? 0.f : 1.f / piPred;
		Float invetar = (etar == 0.f) ? 0.f : 1.f / etar;

		Float ratioDirect = (i == 1) ? state[EDIR] : 1.f;
		Float ratioDirectPred = (i == 2) ? statePred[EDIR] : 1.f;

		Float wVC = ratioDirectPred + MisHeuristic(pir2) * statePred[EVC];
		wVC = ratioDirect + MisHeuristic(pir1 * invpiPred) * wVC;
		wVC *= MisHeuristic(invpi * invetar);

		Float wPM = MisHeuristic(pir2 * etar2) + MisHeuristic(pir2) * statePred[EPM];
		wPM = MisHeuristic(pir1 * etar1) + MisHeuristic(pir1 * invpiPred) * wPM;
		wPM *= MisHeuristic(invpi * invetar);

		return (useVC ? wVC : 0.f) + (useVM ? wPM : 0.f);
	}
	void updateMisHelper(int i, const Path &path, MisState &state, const Scene* scene,
		Float gatherRadius, size_t numLightPath, bool useVC, bool useVM, ETransportMode mode){

		if (i == 0){
			state[EDIR] = 1.f;
		}
		else if (i == 1){
			PathVertex *v0 = path.vertex(1);
			PathVertex *v1 = path.vertex(2);
			EMeasure measure = v0->getAbstractEmitter()->getDirectMeasure();
			Float p0dir = v1->evalPdfDirect(scene, v0, mode, measure == ESolidAngle ? EArea : measure);
			Float p0 = path.vertex(0)->pdf[mode];
			state[EDIR] = MisHeuristic(p0dir / p0);
		}

		if (i == 0)
			state[EVC] = 0.f;
		else{
			PathVertex *v = path.vertex(i + 1);
			PathVertex *vPred = path.vertex(i);
			PathVertex *vPred2 = path.vertex(i - 1);
			Float piPred = vPred2->pdf[mode];
			Float pir2 = vPred->pdf[1 - mode];
			state[EVC] = misEVC(i, state, piPred, pir2, vPred, vPred2);
		}

		if (i < 3)
			state[EPM] = 0.f;
		else{
			PathVertex *v = path.vertex(i + 1);
			PathVertex *vPred = path.vertex(i);
			PathVertex *vPred2 = path.vertex(i - 1);
			PathVertex *vPred3 = path.vertex(i - 2);
			Float piPred = vPred2->pdf[mode];
			Float pir2 = vPred->pdf[1 - mode];
			state[EPM] = misEPM(i, state, piPred, pir2, vPred, vPred3, gatherRadius, numLightPath, mode);
		}
	}
	void initializeMisHelper(const Path &path, MisState *states, const Scene* scene,
		Float gatherRadius, size_t nLightPaths, bool useVC, bool useVM,
		ETransportMode mode){
		for (int i = 0; i < (int)path.vertexCount() - 1; ++i) {
			if (i > 0) states[i] = states[i - 1];
			updateMisHelper(i, path, states[i], scene, gatherRadius, nLightPaths, useVC, useVM, mode);
		}
	}
	Float miWeightVC(const Scene *scene, int s, int t, MisState emitterState, MisState sensorState,
		const PathVertex *vsPred3, const PathVertex *vsPred2, const PathVertex *vsPred, const PathVertex *vs,
		const PathVertex *vt, const PathVertex *vtPred, const PathVertex *vtPred2, const PathVertex *vtPred3,
		Float gatherRadius, size_t numLightPath, bool useVC, bool useVM){

		Float wLight = 0.f, wCamera = 0.f;
		if (s == 0){
			const PositionSamplingRecord &pRec = vt->getPositionSamplingRecord();
			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			EMeasure measure = vt->getAbstractEmitter()->getDirectMeasure();
			Float ptrdir = vtPred->evalPdfDirect(scene, vt, EImportance, measure == ESolidAngle ? EArea : measure);
			Float ptr = evalPdf(vs, scene, NULL, vt, EImportance, measure == ESolidAngle ? EArea : measure);

			Float pt = vtPred->pdf[ERadiance];
			Float ptPred = vtPred2->pdf[ERadiance];
			Float ptr1 = evalPdf(vt, scene, vs, vtPred, EImportance, EArea);
			Float invpt = (pt == 0.f) ? 0.f : 1.f / pt;
			if (useVC)
				wCamera += MisHeuristic(ptrdir * invpt) + MisHeuristic(ptr * ptr1 * invpt) * sensorState[EVC];
			if (useVM)
				wCamera += misWeightVC_pm(t - 1, s - 1, sensorState,
				pt, ptPred, ptr, ptr1,
				vs, vt, vtPred, vtPred2,
				gatherRadius, numLightPath, ERadiance);
		}
		else{
			// MIS weight for light sub path
			EMeasure vsMeasure = EArea;
			Float ratioEmitterDirect = 1.f;
			Float ps = vsPred->pdf[EImportance];
			Float psr = evalPdf(vt, scene, vtPred, vs, ERadiance, EArea);
			if (s == 1){
				const PositionSamplingRecord &pRec = vs->getPositionSamplingRecord();
				const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
				EMeasure measure = emitter->getDirectMeasure();
				Float psdir = vt->evalPdfDirect(scene, vs, EImportance, measure == ESolidAngle ? EArea : measure);
				ratioEmitterDirect = ps / psdir;
				if (!emitter->needsDirectionSample()) // directional light
					vsMeasure = EDiscrete;
			}
			Float psr1 = evalPdf(vs, scene, vt, vsPred, ERadiance, EArea);
			Float psPred = (vsPred2 != NULL) ? vsPred2->pdf[EImportance] : 0.f;
			if (useVC)
				wLight += MisHeuristic(ratioEmitterDirect) * misWeightVC_vc(s - 1, emitterState,
				ps, psr, psr1, vs, vsPred);
			if (useVM)
				wLight += MisHeuristic(ratioEmitterDirect) * misWeightVC_pm(s - 1, t - 1, emitterState,
				ps, psPred, psr, psr1,
				vt, vs, vsPred, vsPred2,
				gatherRadius, numLightPath, EImportance);


			// MIS weight for camera sub path
			if (t == 1){
				wCamera = 0.f;
			}
			else{
				Float ptr = evalPdf(vs, scene, vsPred, vt, EImportance, vsMeasure);
				Float ptr1 = evalPdf(vt, scene, vs, vtPred, EImportance, EArea);
				Float pt = vtPred->pdf[ERadiance];
				//Float ptr2 = (t > 3) ? vtPred->evalPdf(scene, vt, vtPred2, EImportance, EArea) : 0.f;
				Float ptPred = (vtPred2 != NULL) ? vtPred2->pdf[ERadiance] : 0.f;
				if (useVC)
					wCamera += MisHeuristic(ratioEmitterDirect) *										// divide direct sampling prob. correction
					misWeightVC_vc(t - 1, sensorState,
					pt, ptr, ptr1, vt, vtPred);
				if (useVM)
					wCamera += MisHeuristic(ratioEmitterDirect) * misWeightVC_pm(t - 1, s - 1, sensorState,
					pt, ptPred, ptr, ptr1,
					vs, vt, vtPred, vtPred2,
					gatherRadius, numLightPath, ERadiance);
			}
		}

		Float miWeight = 1.f / (wLight + wCamera + 1.f);

		return miWeight;
	}

	Float miWeightVM(const Scene *scene, int s, int t,
		const PathVertex *vsPred2, const PathVertex *vsPred, const PathVertex *vs,
		const PathVertex *vt, const PathVertex *vtPred, const PathVertex *vtPred2,
		MisState emitterState, MisState sensorState,
		bool useVC, bool cameraDirConnection,
		Float gatherRadius = 0.f, size_t numLightPath = 0){
		BDAssert(false); //discarded
	}


	Float miWeightVM(const Scene *scene, int s, int t,
		MisState emitterState, MisState sensorState, MisState emitterStatePred, MisState sensorStatePred,
		const PathVertex *vsPred3, const PathVertex *vsPred2, const PathVertex *vsPred, const PathVertex *vs,
		const PathVertex *vt, const PathVertex *vtPred, const PathVertex *vtPred2, const PathVertex *vtPred3,
		bool cameraDirConnection, Float gatherRadius, size_t numLightPath, bool useVC, bool useVM){

		Float wLight = 0.f, wCamera = 0.f;
		if (cameraDirConnection){
			Float pt = evalPdf(vtPred, scene, vtPred2, vs, ERadiance, EArea);
			Float ps = vsPred->pdf[EImportance];

			Float psPred = vsPred2->pdf[EImportance];
			Float psr1 = evalPdf(vs, scene, vtPred, vsPred, ERadiance, EArea);
			wLight += misWeightPM(s - 1, t - 1, emitterState,
				ps, psPred, pt, psr1,
				vtPred, vs, vsPred, vsPred2,
				gatherRadius, numLightPath, EImportance, useVC, useVM);

			Float ptr1 = evalPdf(vs, scene, vsPred, vtPred, EImportance, EArea);
			Float ptr2 = evalPdf(vtPred, scene, vs, vtPred2, EImportance, EArea);
			Float ptPred = vtPred2->pdf[ERadiance];
			Float ptPred2 = (vtPred3 != NULL) ? vtPred3->pdf[ERadiance] : 0.f;
			wCamera += misWeightPM_pred(t - 1, s - 1, sensorState, sensorStatePred,
				pt, ptPred, ptPred2, ps, ptr1, ptr2,
				vsPred, vs, vtPred, vtPred2, vtPred3,
				gatherRadius, numLightPath, ERadiance, useVC, useVM);
			// 		wCamera += misWeightPM(t - 1, s - 1, sensorState,
			// 			pt, ptPred, ps, ptr1,
			// 			vsPred, vs, vtPred, vtPred2,
			// 			gatherRadius, numLightPath, ERadiance, useVC, useVM);
		}
		else{
			Float ps = evalPdf(vsPred, scene, vsPred2, vt, EImportance, EArea);
			Float pt = vtPred->pdf[ERadiance];

			Float psr1 = evalPdf(vt, scene, vtPred, vsPred, ERadiance, EArea);
			Float psr2 = evalPdf(vsPred, scene, vt, vsPred2, ERadiance, EArea);
			Float psPred = vsPred2->pdf[EImportance];
			Float psPred2 = (vsPred3 != NULL) ? vsPred3->pdf[EImportance] : 0.f;
			wLight += misWeightPM_pred(s - 1, t - 1, emitterState, emitterStatePred,
				ps, psPred, psPred2, pt, psr1, psr2,
				vtPred, vt, vsPred, vsPred2, vsPred3,
				gatherRadius, numLightPath, EImportance, useVC, useVM);
			// 		wLight += misWeightPM(s - 1, t - 1, emitterState,
			// 			ps, psPred, pt, psr1,
			// 			vtPred, vt, vsPred, vsPred2,
			// 			gatherRadius, numLightPath, EImportance, useVC, useVM);

			Float ptr1 = evalPdf(vt, scene, vsPred, vtPred, EImportance, EArea);
			Float ptr2 = evalPdf(vtPred, scene, vt, vtPred2, EImportance, EArea);
			Float ptPred = vtPred2->pdf[ERadiance];
			wCamera += misWeightPM(t - 1, s - 1, sensorState,
				pt, ptPred, ps, ptr1,
				vsPred, vt, vtPred, vtPred2,
				gatherRadius, numLightPath, ERadiance, useVC, useVM);
		}

		Float miWeight = 1.f / (1.f + wLight + wCamera);
		if (!_finite(miWeight) || _isnan(miWeight)){
			SLog(EWarn, "Invalid MIS weight %f in miWeightVM", miWeight);
		}
		return miWeight;
	}

	ref<WorkProcessor> clone() const {
		return new GuidedUPMRenderer(m_config);
	}

	MTS_DECLARE_CLASS()
private:
	GuidedUPMConfiguration m_config;
	ref<Scene> m_scene;
	ref<Sensor> m_sensor;
	ref<Film> m_film;
	ref<PathSampler> m_pathSampler;
	ref<Sampler> m_sampler;
	ref<GuidingSamplers> m_guidingSampler;
};

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

GuidedUPMProcess::GuidedUPMProcess(const RenderJob *parent, RenderQueue *queue,
	const GuidedUPMConfiguration &conf) : m_job(parent), m_queue(queue),
		m_config(conf), m_progress(NULL) {
	m_timeoutTimer = new Timer();
	m_refreshTimer = new Timer();
	m_resultMutex = new Mutex();
	m_resultCounter = 0;
	m_workCounter = 0;
	m_refreshTimeout = 1;
}

ref<WorkProcessor> GuidedUPMProcess::createWorkProcessor() const {
	return new GuidedUPMRenderer(m_config);
}

void GuidedUPMProcess::develop() {
	LockGuard lock(m_resultMutex);
	size_t pixelCount = m_result->getImageBlock()->getBitmap()->getPixelCount();
	const Spectrum *accum = (Spectrum *)(m_result->getImageBlock()->getBitmap()->getData());
	Spectrum *target = (Spectrum *) m_developBuffer->getData();
 	Float invFactor = 1.f / Float(m_result->getSampleCount());
	for (size_t i=0; i<pixelCount; ++i) {
		target[i] = accum[i] * invFactor;
	}	
	m_film->setBitmap(m_developBuffer);
	m_refreshTimer->reset();
	m_queue->signalRefresh(m_job);
}

void GuidedUPMProcess::processResult(const WorkResult *workResult, bool cancelled) {
	const UPMWorkResult *wr = static_cast<const UPMWorkResult *>(workResult);
	LockGuard lock(m_resultMutex);
	if (m_resultCounter == 0) m_result->clear();	
	m_result->put(wr);
	m_progress->update(++m_resultCounter);
	m_refreshTimeout = std::min(2000U, m_refreshTimeout * 2);

	/* Re-develop the entire image every two seconds if partial results are
	   visible (e.g. in a graphical user interface). */
	if (m_job->isInteractive() && m_refreshTimer->getMilliseconds() > m_refreshTimeout)
		develop();
}

ParallelProcess::EStatus GuidedUPMProcess::generateWork(WorkUnit *unit, int worker) {
	int64_t timeout = 0;
	if (m_config.timeout > 0) {
		timeout = static_cast<int64_t>(m_config.timeout) - static_cast<int64_t>(m_timeoutTimer->getSeconds());
	}

	if (m_workCounter >= m_config.workUnits || timeout < 0)
		return EFailure;

	SeedWorkUnit *workUnit = static_cast<SeedWorkUnit *>(unit);
	workUnit->setID(m_workCounter++);
	workUnit->setTotalWorkNum(m_config.workUnits);
	workUnit->setTimeout(timeout);
	return ESuccess;
}

void GuidedUPMProcess::bindResource(const std::string &name, int id) {
	ParallelProcess::bindResource(name, id);
	if (name == "sensor") {
		m_film = static_cast<Sensor *>(Scheduler::getInstance()->getResource(id))->getFilm();
		if (m_progress)
			delete m_progress;
		m_progress = new ProgressReporter("Rendering", m_config.workUnits, m_job);
		m_result = new UPMWorkResult(m_film->getCropSize().x, m_film->getCropSize().y, m_config.maxDepth, NULL, true);
		m_result->clear();		
		m_developBuffer = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, m_film->getCropSize());
	}
}

MTS_IMPLEMENT_CLASS_S(GuidedUPMRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(GuidedUPMProcess, false, ParallelProcess)
MTS_IMPLEMENT_CLASS(SeedWorkUnit, false, WorkUnit)
MTS_IMPLEMENT_CLASS(UPMWorkResult, false, WorkResult)
MTS_NAMESPACE_END
