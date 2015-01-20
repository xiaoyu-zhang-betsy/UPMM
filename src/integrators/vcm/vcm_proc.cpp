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

#include <mitsuba/core/sfcurve.h>
#include <mitsuba/core/statistics.h>
#include <mitsuba/bidir/util.h>
#include <mitsuba/bidir/path.h>
#include "vcm_proc.h"


MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                         Worker implementation                        */
/* ==================================================================== */

StatsCounter largeStepRatio("Primary sample space MLT",
	"Accepted large steps", EPercentage);
StatsCounter smallStepRatio("Primary sample space MLT",
	"Accepted small steps", EPercentage);
StatsCounter acceptanceRate("Primary sample space MLT",
	"Overall acceptance rate", EPercentage);
StatsCounter forcedAcceptance("Primary sample space MLT",
	"Number of forced acceptances");

/*
*	Misc for VCM
*/
inline Float MisHeuristic(Float pdf) {
	return pdf * pdf;
}
enum EMisTechV {
	EVCMV = 0,
	EVCV = 1,
	EVMV = 2,
	EMisTechsV = 3
};
struct MisStateV{
	float state[EMisTechs];
	MisStateV(){
		state[EVCMV] = 0.f;
		state[EVCV] = 0.f;
		state[EVMV] = 0.f;
	}
	float& operator[](EMisTechV tech){
		return state[tech];
	}
};
struct LightVertexV{
	Spectrum importanceWeight;
	Vector wo;
	MisStateV emitterState;
	inline  LightVertexV(const PathVertex *vs, const PathVertex *vsPred,
		const MisStateV &state, Spectrum _wgt) :
		emitterState(state), importanceWeight(_wgt){
		wo = normalize(vsPred->getPosition() - vs->getPosition());
	}
};
struct LightVertexExtV{
	int depth;
	const Shape *shape;
	Point3 position;
	Frame shFrame;
	Frame geoFrame;
	uint8_t measure;
	uint8_t type;

	bool hasVsPred;
	Point3 posPred;
	Float pdfPred;
	uint8_t typePred;
	uint8_t measPred;

	inline LightVertexExtV(const PathVertex *vs, const PathVertex *vsPred, int _depth){
		const Intersection &its = vs->getIntersection();
		depth = _depth;
		shape = its.shape;
		position = its.p;
		shFrame = its.shFrame;
		geoFrame = its.geoFrame;
		measure = vs->measure;
		type = vs->type;
		if (vsPred != NULL){
			hasVsPred = true;
			posPred = vsPred->getPosition();
			pdfPred = vsPred->pdf[EImportance];
			typePred = vsPred->type;
			measPred = vsPred->measure;
		}
		else{
			hasVsPred = false;
		}
	}
};
struct LightPathNodeDataV{
	int depth;
	size_t vertexIndex;
};
struct LightPathNodeV : public SimpleKDNode < Point, LightPathNodeDataV > {

	inline LightPathNodeV(){}

	inline  LightPathNodeV(const Point3 p, size_t vertexIndex, int depth){
		position = p;
		data.vertexIndex = vertexIndex;
		data.depth = depth;
	}
	inline LightPathNodeV(Stream *stream) {
		// TODO
	}
	void serialize(Stream *stream) const {
		// TODO
	}
};
typedef PointKDTree<LightPathNodeV>		LightPathTreeV;
typedef LightPathTreeV::IndexType     IndexTypeV;
typedef LightPathTreeV::SearchResult SearchResultV;

struct VertexMergingQuery {
	VertexMergingQuery(const Scene *_scene,
		const PathVertex *_vt, const PathVertex *_vtPred, const Vector _wi, const Vector _wiPred,
		const Spectrum &_radianceWeight, const std::vector<LightVertexV> &lightVertices,
		int _t, int _maxDepth, Float _misVcWeightFactor, Float _vmNormalization,
		const MisStateV &_sensorState)
		: scene(_scene), radianceWeight(_radianceWeight), t(_t), maxDepth(_maxDepth),
		vt(_vt), vtPred(_vtPred), wiPred(_wiPred), wi(_wi),
		misVcWeightFactor(_misVcWeightFactor), vmNormalization(_vmNormalization),
		sensorState(_sensorState), m_lightVertices(lightVertices),
		result(0.f){}

	inline void operator()(const LightPathNodeV &path) {
		int s = path.data.depth;
		if (maxDepth != -1 && s + t > maxDepth + 2) return;

		size_t vertexIndex = path.data.vertexIndex;
		LightVertexV v = m_lightVertices[vertexIndex];
		Vector wo = v.wo;
		MisStateV emitterState = v.emitterState;
		Spectrum bsdfFactor = vt->eval(scene, wi, wo, ERadiance);
		Spectrum val = v.importanceWeight * radianceWeight * bsdfFactor;

		if (val.isZero()) return;

		Float weightExt = 0.f;
		Float psr2_w = vt->evalPdf(scene, wi, wo, ERadiance, ESolidAngle);
		Float ptr2_w = vt->evalPdf(scene, wo, wi, EImportance, ESolidAngle);
		Float wLight = emitterState[EVCMV] * MisHeuristic(misVcWeightFactor) + MisHeuristic(psr2_w) * emitterState[EVMV];
		Float wCamera = sensorState[EVCMV] * MisHeuristic(misVcWeightFactor) + MisHeuristic(ptr2_w) * sensorState[EVMV];
		weightExt = 1.f / (1.f + wLight + wCamera);

		result += val * weightExt * vmNormalization;
	}
	int t, maxDepth;
	Float misVcWeightFactor, vmNormalization;
	MisStateV sensorState;
	Vector wi, wiPred;
	const Scene *scene;
	const PathVertex *vt, *vtPred;
	const Spectrum &radianceWeight;
	const std::vector<LightVertexV> &m_lightVertices;
	Spectrum result;
};


class VCMRenderer : public WorkProcessor {
public:
	VCMRenderer(const VCMConfiguration &conf)
		: m_config(conf) {
	}

	VCMRenderer(Stream *stream, InstanceManager *manager)
		: WorkProcessor(stream, manager) {
		m_config = VCMConfiguration(stream);
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		m_config.serialize(stream);
	}

	ref<WorkUnit> createWorkUnit() const {
		return new SeedWorkUnit();
	}

	ref<WorkResult> createWorkResult() const {
		return new UPMWorkResult(m_film->getCropSize().x, m_film->getCropSize().y, m_config.maxDepth, m_film->getReconstructionFilter());
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
	}

	void updateMisHelper(int i, const Path &path, MisStateV &state, const Scene* scene,
		size_t nLightPaths, Float misVcWeightFactor, Float misVmWeightFactor, ETransportMode mode){
		if (i == 0){
			state[EVCMV] = 0.f;
			state[EVCV] = 0.f;
			state[EVMV] = 0.f;
		}
		else if (i == 1){
			state[EVCMV] = 0.f;
			state[EVCV] = 0.f;
			state[EVMV] = 0.f;
			PathVertex *v = path.vertex(1);
			PathEdge *e = path.edge(1);
			PathVertex *vNext = path.vertex(2);
			EMeasure measure = v->getAbstractEmitter()->getDirectMeasure();
			Float pconnect = vNext->evalPdfDirect(scene, v, mode, measure == ESolidAngle ? EArea : measure);
			Float ptrace = path.vertex(0)->pdf[mode] * path.edge(0)->pdf[mode];
			Float p1 = v->pdf[mode] * e->pdf[mode];
			if (mode == ERadiance){
				// 			if (measure == EDiscrete)
				// 				p1 *= (e->length * e->length) / (vNext->isOnSurface() ? std::abs(dot(e->d, vNext->getGeometricNormal())) : 1.f);
				state[EVCMV] = MisHeuristic(pconnect / (ptrace * p1));
				state[EVCV] = 0.f;
				state[EVMV] = 0.f;
			}
			else{
				state[EVCMV] = MisHeuristic(pconnect / (ptrace * p1));
				if (measure != EDiscrete){
					Float geoTerm0 = v->isOnSurface() ? std::abs(dot(e->d, v->getGeometricNormal()) / (e->length * e->length)) : 1;
					state[EVCV] = MisHeuristic(geoTerm0 / (ptrace * p1));
				}
				else{
					state[EVCV] = 0.f;
				}
				state[EVMV] = state[EVCV] * MisHeuristic(misVcWeightFactor);
			}
		}
		else{
			PathVertex *v = path.vertex(i);
			PathEdge   *e = path.edge(i);
			if (v->measure == EDiscrete){
				PathVertex *vNext = path.vertex(i + 1);
				state[EVCMV] = 0.f;
				state[EVCV] *= MisHeuristic(std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / std::abs(vNext->isOnSurface() ? dot(e->d, vNext->getGeometricNormal()) : 1));
				state[EVMV] *= MisHeuristic(std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / std::abs(vNext->isOnSurface() ? dot(e->d, vNext->getGeometricNormal()) : 1));
			}
			else{
				BDAssert(v->measure == EArea);
				PathVertex *vPred = path.vertex(i - 1);
				PathEdge * ePred = path.edge(i - 1);
				Float giIn = std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / (e->length * e->length);
				Float giInPred = std::abs(vPred->isOnSurface() ? dot(ePred->d, vPred->getGeometricNormal()) : 1) / (ePred->length * ePred->length);
				Float pi = v->pdf[mode] * e->pdf[mode];
				Float pir2_w = v->pdf[1 - mode] * ePred->pdf[1 - mode] / giInPred;
				state[EVCV] = MisHeuristic(giIn / pi) * (state[EVCMV] + MisHeuristic(misVmWeightFactor) + MisHeuristic(pir2_w) * state[EVCV]);
				state[EVMV] = MisHeuristic(giIn / pi) * (1.f + state[EVCMV] * MisHeuristic(misVcWeightFactor) + MisHeuristic(pir2_w) * state[EVMV]);
				state[EVCMV] = MisHeuristic(1 / pi);
			}
		}
	}
	void initializeMisHelper(const Path &path, MisStateV *states, const Scene* scene,
		size_t nLightPaths, Float misVcWeightFactor, Float misVmWeightFactor, ETransportMode mode){
		for (int i = 0; i < (int)path.vertexCount() - 1; ++i) {
			if (i > 0) states[i] = states[i - 1];
			updateMisHelper(i, path, states[i], scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, mode);
		}
	}

	Float miWeightVC(const Scene *scene,
		const PathVertex *vsPred, const PathVertex *vs,
		const PathVertex *vt, const PathVertex *vtPred,
		int s, int t, bool sampleDirect,
		float emitterdVCM, float emitterdVC,
		float sensordVCM, float sensordVC,
		Float misVmWeightFactor, size_t nLightPaths){
		Float weight = 0.0f;

		if (s == 0){
			const PositionSamplingRecord &pRec = vt->getPositionSamplingRecord();
			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			EMeasure measure = vt->getAbstractEmitter()->getDirectMeasure();
			Float ptrdir = vtPred->evalPdfDirect(scene, vt, EImportance, measure == ESolidAngle ? EArea : measure);
			Float ptr = vs->evalPdf(scene, NULL, vt, EImportance, measure == ESolidAngle ? EArea : measure);
			Float ptr2_w = vt->evalPdf(scene, vs, vtPred, EImportance, ESolidAngle);
			//Float wCamera = (t == 1) ? 0.f : ratioEmitterDirect * ptr1 * (misVmWeightFactor + sensordVCM + ptr2_w * sensordVC);
			Float wCamera = (t == 1) ? 0.f : MisHeuristic(ptrdir) * sensordVCM + MisHeuristic(ptr * ptr2_w) * sensordVC;
			weight = 1.f / (1.f + wCamera);
		}
		else if (s > 1 && t > 1){
			Float psr2_w = vs->evalPdf(scene, vt, vsPred, ERadiance, ESolidAngle);
			Float ptr2_w = vt->evalPdf(scene, vs, vtPred, EImportance, ESolidAngle);
			Float psr1 = vt->evalPdf(scene, vtPred, vs, ERadiance, EArea);
			Float ptr1 = vs->evalPdf(scene, vsPred, vt, EImportance, EArea);
			Float wLight = MisHeuristic(psr1) * (MisHeuristic(misVmWeightFactor) + emitterdVCM + MisHeuristic(psr2_w) * emitterdVC);
			Float wCamera = MisHeuristic(ptr1) * (MisHeuristic(misVmWeightFactor) + sensordVCM + MisHeuristic(ptr2_w) * sensordVC);
			weight = 1.f / (1.f + wLight + wCamera);
		}
		else if (t == 1 && s > 1){
			EMeasure measure = vt->getAbstractEmitter()->getDirectMeasure();
			Float pconnect = vs->evalPdfDirect(scene, vt, ERadiance, measure == ESolidAngle ? EArea : measure);
			Float ptrace = vtPred->pdf[ERadiance];
			Float psr2_w = vs->evalPdf(scene, vt, vsPred, ERadiance, ESolidAngle);
			Float psr1 = vt->evalPdf(scene, vtPred, vs, ERadiance, EArea);

			Float wLight = MisHeuristic(ptrace * psr1 / pconnect) * (MisHeuristic(misVmWeightFactor) + emitterdVCM + MisHeuristic(psr2_w) * emitterdVC);
			weight = 1.f / (1.f + wLight);
		}
		else if (s == 1){
			EMeasure vsMeasure = EArea;
			Float ratioEmitterDirect = 1.f;
			Float ptrace = vsPred->pdf[EImportance];
			const PositionSamplingRecord &pRec = vs->getPositionSamplingRecord();
			const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
			if (s == 1){				
				EMeasure measure = emitter->getDirectMeasure();
				Float psdir = vt->evalPdfDirect(scene, vs, EImportance, measure == ESolidAngle ? EArea : measure);
				ratioEmitterDirect = ptrace / psdir;
				if (!emitter->needsDirectionSample()) // directional light
					vsMeasure = EDiscrete;
			}

			Float ptr2_w = vt->evalPdf(scene, vs, vtPred, EImportance, ESolidAngle);
			Float psr1 = (!emitter->needsDirectionSample() || !emitter->needsPositionSample()) ? 0.f : vt->evalPdf(scene, vtPred, vs, ERadiance, EArea);
			Float ptr1 = vs->evalPdf(scene, vsPred, vt, EImportance, vsMeasure);
			Float wLight = MisHeuristic(ratioEmitterDirect * (psr1 / ptrace));
			Float wCamera = (t == 1) ? 0.f : MisHeuristic(ratioEmitterDirect * ptr1) * (MisHeuristic(misVmWeightFactor) + sensordVCM + MisHeuristic(ptr2_w) * sensordVC);
			weight = 1.f / (1.f + wLight + wCamera);
		}
		else{
			SLog(EError, "Not implemented miweightVC for s = %d, t = %d", s, t);
		}
		return weight;
	}

	void gatherLightPaths(ref<PathSampler> pathSampler,
		const bool useVC, const bool useVM,
		const float gatherRadius, const int nsample, 
		UPMWorkResult *wr, ImageBlock *batres){
		const Sensor *sensor = m_scene->getSensor();
		m_lightPathTree.clear();
		m_lightVertices.clear();
		m_lightVerticesExt.clear();
		m_lightPathEnds.clear();

		PathVertex vtPred, vt;
		PathEdge vtEdge, connectionEdge;

		Float time = sensor->getShutterOpen();
		Vector2i filmSize = sensor->getFilm()->getSize();
		m_lightPathNum = nsample;
		Float etaVCM = (M_PI * gatherRadius * gatherRadius) * m_lightPathNum;
		Float invLightPathNum = 1.f / m_lightPathNum;
		Float misVmWeightFactor = useVM ? etaVCM : 0.f;
		Float misVcWeightFactor = useVC ? 1.f / etaVCM : 0.f;
		for (size_t k = 0; k < m_lightPathNum; k++){
			/* Initialize the path endpoints */
			pathSampler->m_emitterSubpath.initialize(m_scene, time, EImportance, pathSampler->m_pool);

			/* Perform random walks from the emitter side */
			pathSampler->m_emitterSubpath.randomWalk(m_scene, m_sampler, pathSampler->m_emitterDepth,
				pathSampler->m_rrDepth, EImportance, pathSampler->m_pool);

			// emitter states
			MisStateV emitterState, sensorState;
			Spectrum importanceWeight = Spectrum(1.0f);
			for (int s = 1; s < (int)pathSampler->m_emitterSubpath.vertexCount(); ++s) {
				PathVertex
					*vsPred2 = pathSampler->m_emitterSubpath.vertexOrNull(s - 2),
					*vsPred = pathSampler->m_emitterSubpath.vertex(s - 1),
					*vs = pathSampler->m_emitterSubpath.vertex(s);
				PathEdge
					*esPred = pathSampler->m_emitterSubpath.edgeOrNull(s - 2),
					*es = pathSampler->m_emitterSubpath.edgeOrNull(s - 1);

				// Compute the throughput of emitter subpath till this vertex
				importanceWeight *= vsPred->weight[EImportance] * vsPred->rrWeight * es->weight[EImportance];

				// update mis helper and path type
				updateMisHelper(s - 1, pathSampler->m_emitterSubpath, emitterState, m_scene, m_lightPathNum, misVcWeightFactor, misVmWeightFactor, EImportance);

				// store light paths												
				if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon /* don't save backfaced photons */){
					LightVertexV lvertex(vs, vsPred, emitterState, importanceWeight);
					m_lightVertices.push_back(lvertex);
					LightVertexExtV lvertexExt = LightVertexExtV(vs, vsPred, s);
					m_lightVerticesExt.push_back(lvertexExt);
					LightPathNodeV lnode(vs->getPosition(), m_lightVertices.size() - 1, s);
					m_lightPathTree.push_back(lnode);
				}

				// connect to camera
				if (vs->measure != EDiscrete && wr != NULL && useVC){
					Spectrum value = importanceWeight * vs->sampleDirect(m_scene, pathSampler->m_directSampler,
						&vtPred, &vtEdge, &vt, ERadiance);
					if (value.isZero()) continue;

					value *= vs->eval(m_scene, vsPred, &vt, EImportance);
					vs->measure = EArea;
					int interactions = pathSampler->m_maxDepth - s;
					if (value.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, &vt, &vtEdge, interactions))
						continue;
					value *= connectionEdge.evalCached(vs, &vt, PathEdge::ETransmittance | PathEdge::ECosineImp);
					Float weight = miWeightVC(m_scene, vsPred, vs, &vt, &vtPred,
						s, 1, pathSampler->m_sampleDirect,
						emitterState[EVCMV], emitterState[EVCV],
						sensorState[EVCMV], sensorState[EVCV],
						misVmWeightFactor, m_lightPathNum);
					value *= weight;
					if (value.isZero()) continue;

					Point2 samplePos(0.0f);
					vt.getSamplePosition(vs, samplePos);

#if UPM_DEBUG == 1
					wr->putDebugSampleVC(samplePos, value);
#endif					
					wr->putSample(samplePos, &value[0]);
					if (batres != NULL)
						batres->put(samplePos, &value[0]);
				}
			}
			m_lightPathEnds.push_back(m_lightVertices.size());
		}
		// build kdtree
		m_lightPathTree.build(true);
		/* Release any used edges and vertices back to the memory pool */
		pathSampler->m_emitterSubpath.release(pathSampler->m_pool);
	}
	void sampleCameraPath(ref<PathSampler> pathSampler, UPMWorkResult *wr, 
		const bool useVC, const bool useVM,
		const float gatherRadius, const Point2i &offset, const size_t cameraPathIndex, SplatList &list) {
		list.clear();

		const Sensor *sensor = m_scene->getSensor();
		size_t nLightPaths = m_lightPathNum;
		const float etaVCM = (M_PI * gatherRadius * gatherRadius) * nLightPaths;
		Float vmNormalization = 1.f / etaVCM;
		Float misVmWeightFactor = useVM ? etaVCM : 0.f;
		Float misVcWeightFactor = useVC ? 1.f / etaVCM : 0.f;
		Vector2i filmSize = sensor->getFilm()->getSize();

		switch (pathSampler->m_technique) {
		case PathSampler::EBidirectional: {
			/* Uniformly sample a scene time */
			Float time = sensor->getShutterOpen();
			if (sensor->needsTimeSample())
				time = sensor->sampleTime(pathSampler->m_sensorSampler->next1D());

			/* Initialize the path endpoints */
			pathSampler->m_sensorSubpath.initialize(m_scene, time, ERadiance, pathSampler->m_pool);

			if (offset == Point2i(-1))
				pathSampler->m_sensorSubpath.randomWalk(m_scene, pathSampler->m_sensorSampler,
					pathSampler->m_sensorDepth, pathSampler->m_rrDepth, ERadiance, pathSampler->m_pool);
			else
				pathSampler->m_sensorSubpath.randomWalkFromPixel(m_scene, pathSampler->m_sensorSampler,
					pathSampler->m_sensorDepth, offset, pathSampler->m_rrDepth, pathSampler->m_pool);

			/* Compute the combined weights along the two subpaths */
			Spectrum *radianceWeights = (Spectrum *)alloca(pathSampler->m_sensorSubpath.vertexCount()  * sizeof(Spectrum));

			radianceWeights[0] = Spectrum(1.0f);

			for (size_t i = 1; i<pathSampler->m_sensorSubpath.vertexCount(); ++i)
				radianceWeights[i] = radianceWeights[i - 1] *
				pathSampler->m_sensorSubpath.vertex(i - 1)->weight[ERadiance] *
				pathSampler->m_sensorSubpath.vertex(i - 1)->rrWeight *
				pathSampler->m_sensorSubpath.edge(i - 1)->weight[ERadiance];

			Point2 initialSamplePos;
			if (pathSampler->m_sensorSubpath.vertexCount() > 2) {
				Point2 samplePos(0.0f);
				pathSampler->m_sensorSubpath.vertex(1)->getSamplePosition(pathSampler->m_sensorSubpath.vertex(2), samplePos);
				list.append(samplePos, Spectrum(0.0f));
				initialSamplePos = samplePos;
			}

			// initialize of MIS helper
			// VCM tech report Eq. 31 - 33		
			MisStateV *sensorStates = (MisStateV *)alloca(pathSampler->m_sensorSubpath.vertexCount() * sizeof(MisStateV));
			initializeMisHelper(pathSampler->m_sensorSubpath, sensorStates, m_scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, ERadiance);

			PathVertex tempEndpoint, tempSample;
			PathEdge tempEdge, connectionEdge;
			Point2 samplePos(0.0f);

			// massive vertex merging
			if (useVM){
				int minT = 2, maxT = (int)pathSampler->m_sensorSubpath.vertexCount() - 1;
				if (pathSampler->m_maxDepth != -1)
					maxT = std::min(maxT, pathSampler->m_maxDepth - 1);

				for (int t = minT; t <= maxT; ++t) {
					PathVertex
						*vtPred2 = pathSampler->m_sensorSubpath.vertexOrNull(t - 2),
						*vtPred = pathSampler->m_sensorSubpath.vertexOrNull(t - 1),
						*vt = pathSampler->m_sensorSubpath.vertex(t);
					PathEdge
						*vtEdge = pathSampler->m_sensorSubpath.edgeOrNull(t - 1);

					if (!vt->isDegenerate()){
						BDAssert(vt->type == PathVertex::ESurfaceInteraction);

						Vector wi = normalize(vtPred->getPosition() - vt->getPosition());
						Vector wiPred = (t == 2) ? Vector(-1.f, -1.f, -1.f) : normalize(vtPred2->getPosition() - vtPred->getPosition());
						VertexMergingQuery query(m_scene, vt, vtPred, wi, wiPred, radianceWeights[t], m_lightVertices,
							t, pathSampler->m_maxDepth, misVcWeightFactor, vmNormalization, sensorStates[t - 1]);

						m_lightPathTree.executeQuery(vt->getPosition(), gatherRadius, query);

#if UPM_DEBUG == 1
						wr->putDebugSampleVM(initialSamplePos, query.result);
#endif

						if (!query.result.isZero())
							list.accum(0, query.result);
					}
				}
			}

			// vertex connection
			if (useVC){
				PathVertex *vs0 = pathSampler->m_pool.allocVertex();
				PathVertex *vsPred0 = pathSampler->m_pool.allocVertex();
				PathEdge *vsEdge0 = pathSampler->m_pool.allocEdge();
				size_t lightPathBegin = (cameraPathIndex == 0) ? 0 : m_lightPathEnds[cameraPathIndex - 1];
				size_t lightPathEnd = m_lightPathEnds[cameraPathIndex];
				for (size_t i = lightPathBegin; i < lightPathEnd + 2; i++){
					int s = lightPathEnd + 1 - i;
					PathVertex* vsPred = NULL, *vs = NULL;
					PathEdge *vsEdge = NULL;
					Spectrum importanceWeight = Spectrum(1.0f);
					MisStateV emitterState;

					if (i < lightPathEnd){
						LightVertexV lvertex = m_lightVertices[i];
						importanceWeight = lvertex.importanceWeight;
						emitterState = lvertex.emitterState;
						LightVertexExtV lvertexExt = m_lightVerticesExt[i];
						s = lvertexExt.depth;
						vs = vs0;
						Intersection &itp = vs->getIntersection();
						itp.p = lvertexExt.position;
						itp.shFrame = lvertexExt.shFrame;
						itp.geoFrame = lvertexExt.geoFrame;
						itp.setShapePointer(lvertexExt.shape);
						vs->measure = lvertexExt.measure;
						vs->type = lvertexExt.type;
						if (lvertexExt.hasVsPred){
							vsPred = vsPred0;
							vsPred->type = lvertexExt.typePred;
							vsPred->measure = lvertexExt.measPred;
							vsPred->setPosition(lvertexExt.posPred);
							vsPred->pdf[EImportance] = lvertexExt.pdfPred;
						}
						else
							vsPred = NULL;
					}

					int minT = 2, maxT = (int)pathSampler->m_sensorSubpath.vertexCount() - 1;
					if (pathSampler->m_maxDepth != -1)
						maxT = std::min(maxT, pathSampler->m_maxDepth + 1 - s);
					for (int t = minT; t <= maxT; ++t) {
						PathVertex
							*vtPred = pathSampler->m_sensorSubpath.vertexOrNull(t - 1),
							*vt = pathSampler->m_sensorSubpath.vertex(t);
						PathEdge
							*vtEdge = pathSampler->m_sensorSubpath.edgeOrNull(t - 1);

						/* Will receive the path weight of the (s, t)-connection */
						Spectrum value = Spectrum(0.0f);

						if (s == 0){
							/* If possible, convert 'vt' into an emitter sample */
							if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
								continue;
							vs = vs0;
							vs->makeEndpoint(m_scene, time, EImportance);
						}
						else if (s == 1){
							value = radianceWeights[t] * vt->sampleDirect(m_scene, pathSampler->m_directSampler,
								&tempEndpoint, &tempEdge, &tempSample, EImportance);
							if (value.isZero())
								continue;
							vs = &tempSample; vsPred = &tempEndpoint; vsEdge = &tempEdge;
						}

						/* Determine the pixel sample position when necessary */
						if (vt->isSensorSample() && !vt->getSamplePosition(vs, samplePos))
							continue;

						RestoreMeasureHelper rmh0(vs), rmh1(vt);

						/* Will be set to true if direct sampling was used */
						bool sampleDirect = false;

						/* Number of edges of the combined subpaths */
						int depth = s + t - 1;

						/* Allowed remaining number of ENull vertices that can
						be bridged via pathConnect (negative=arbitrarily many) */
						int remaining = pathSampler->m_maxDepth - depth;

						// backup original path vertex measure
						uint8_t vsMeasure0 = vs->measure, vtMeasure0 = vt->measure;

						/* Account for the terms of the measurement contribution
						function that are coupled to the connection endpoints */
						if (s == 0){
							/* If possible, convert 'vt' into an emitter sample */
							value = radianceWeights[t] *
								vs->eval(m_scene, NULL, vt, EImportance) *
								vt->eval(m_scene, vtPred, vs, ERadiance);
						}
						else if (pathSampler->m_sampleDirect && ((t == 1 && s > 1) || (s == 1 && t > 1))) {
							if (s == 1) {
								/* Generate a position on an emitter using direct sampling */
								// 							value = radianceWeights[t] * vt->sampleDirect(m_scene, m_directSampler,
								// 								&tempEndpoint, &tempEdge, &tempSample, EImportance);
								// 							if (value.isZero())
								// 								continue;
								// 							vs = &tempSample; vsPred = &tempEndpoint; vsEdge = &tempEdge;
								value *= vt->eval(m_scene, vtPred, vs, ERadiance);
								vt->measure = EArea;
							}
							else {
								/* s==1/t==1 path: use a direct sampling strategy if requested */
								if (vs->isDegenerate())
									continue;
								/* Generate a position on the sensor using direct sampling */
								value = importanceWeight * vs->sampleDirect(m_scene, pathSampler->m_directSampler,
									&tempEndpoint, &tempEdge, &tempSample, ERadiance);
								if (value.isZero())
									continue;
								vt = &tempSample; vtPred = &tempEndpoint; vtEdge = &tempEdge;
								value *= vs->eval(m_scene, vsPred, vt, EImportance);
								vs->measure = EArea;
							}
							sampleDirect = true;
						}
						else {
							/* Can't connect degenerate endpoints */
							if ((vs->isDegenerate()) || vt->isDegenerate())
								continue;
							Intersection &itpt = vs->getIntersection();
							Intersection &itpt2 = vsPred->getIntersection();

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

						if (value.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, vt, vtEdge, interactions))
							continue;

						depth += interactions;

						if (pathSampler->m_excludeDirectIllum && depth <= 2)
							continue;

						/* Account for the terms of the measurement contribution
						function that are coupled to the connection edge */
						if (!sampleDirect)
							value *= connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
						else
							value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
							(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

						MisStateV sensorState = sensorStates[t - 1];
						Float weight = miWeightVC(m_scene, vsPred, vs, vt, vtPred,
							s, t, pathSampler->m_sampleDirect,
							emitterState[EVCMV], emitterState[EVCV],
							sensorState[EVCMV], sensorState[EVCV],
							misVmWeightFactor, nLightPaths);
						value *= weight;

						if (weight == 0.0f) continue;

						if (t < 2) {
							list.append(samplePos, value);
#if UPM_DEBUG == 1
							wr->putDebugSampleVC(samplePos, value);
#endif
						}
						else {
							BDAssert(pathSampler->m_sensorSubpath.vertexCount() > 2);
							list.accum(0, value);
#if UPM_DEBUG == 1
							wr->putDebugSampleVC(initialSamplePos, value);
#endif
						}
					}
				}
				pathSampler->m_pool.release(vs0);
				pathSampler->m_pool.release(vsPred0);
				pathSampler->m_pool.release(vsEdge0);
			}
			else{
				// complement missing path in single VM technique set
				int t = 2, s = 0;
				PathVertex
					*vtPred = pathSampler->m_sensorSubpath.vertexOrNull(t - 1),
					*vt = pathSampler->m_sensorSubpath.vertex(t);
				PathEdge
					*vtEdge = pathSampler->m_sensorSubpath.edgeOrNull(t - 1);
				PathVertex *vs = pathSampler->m_pool.allocVertex();

				/* Will receive the path weight of the (s, t)-connection */
				Spectrum value = Spectrum(0.0f);

				/* If possible, convert 'vt' into an emitter sample */
				if (vt->cast(m_scene, PathVertex::EEmitterSample) && !vt->isDegenerate() && !pathSampler->m_excludeDirectIllum){
					vs->makeEndpoint(m_scene, time, EImportance);

					/* Number of edges of the combined subpaths */
					int depth = s + t - 1;

					/* Allowed remaining number of ENull vertices that can
					be bridged via pathConnect (negative=arbitrarily many) */
					int remaining = pathSampler->m_maxDepth - depth;

					// backup original path vertex measure
					uint8_t vsMeasure0 = vs->measure, vtMeasure0 = vt->measure;

					value = radianceWeights[t] *
						vs->eval(m_scene, NULL, vt, EImportance) *
						vt->eval(m_scene, vtPred, vs, ERadiance);

					/* Attempt to connect the two endpoints, which could result in
					the creation of additional vertices (index-matched boundaries etc.) */
					int interactions = remaining;
					connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, vt, vtEdge, interactions);
					depth += interactions;

					value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
						(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

					MisStateV sensorState = sensorStates[t - 1];
					MisStateV emitterState;
					Float weight = miWeightVC(m_scene, NULL, vs, vt, vtPred,
						s, t, pathSampler->m_sampleDirect,
						emitterState[EVCMV], emitterState[EVCV],
						sensorState[EVCMV], sensorState[EVCV],
						misVmWeightFactor, nLightPaths);
					value *= weight;
					list.accum(0, value);
#if UPM_DEBUG == 1
					wr->putDebugSampleVC(initialSamplePos, value);
#endif
				}
				pathSampler->m_pool.release(vs);
			}

			/* Release any used edges and vertices back to the memory pool */
			pathSampler->m_sensorSubpath.release(pathSampler->m_pool);
		}
			break;

		default:
			Log(EError, "PathSampler::sample(): invalid technique!");
		}
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {				
		UPMWorkResult *result = static_cast<UPMWorkResult *>(workResult);		
		const SeedWorkUnit *wu = static_cast<const SeedWorkUnit *>(workUnit);
		const int workID = wu->getID();
		SplatList *splats = new SplatList();		

		HilbertCurve2D<int> hilbertCurve;
		TVector2<int> filmSize(m_film->getCropSize());
		hilbertCurve.initialize(filmSize);
		uint64_t iteration = workID;
		size_t actualSampleCount = 0;
		ImageBlock *batres = NULL;

#if UPM_DEBUG == 1
		// [UC] for unbiased check		
		Float sepInterval = 60.f;
		size_t numSepSamples = 0;
		size_t numBatch = 0;
		ref<Timer> intervalTimer = new Timer(false);
		if (m_config.enableSeparateDump){
			batres = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize(), m_film->getReconstructionFilter());
			batres->clear();
			intervalTimer->start();
		}

		// [UC] for relative contribution graphing
		Float rcInterval = 600.f;
		Float numProgressiveBatch = 0;
		ref<Timer> intervalTimer2 = new Timer(false);
		if (m_config.enableProgressiveDump){
			intervalTimer2->start();
		}
#endif

		float radius = m_config.initialRadius;
		ref<Timer> timer = new Timer();
		for (size_t j = 0; j < m_config.sampleCount || (wu->getTimeout() > 0 && (int)timer->getMilliseconds() < wu->getTimeout()); j++) {
			if (m_config.initialRadius > 0.0f){
				Float reduceFactor = 1.f / std::pow(Float(iteration + 1), (Float)(0.5 * (1 - 0.75/*radiusAlpha*/)));
				radius = std::max(reduceFactor * m_config.initialRadius, (Float)1e-7);
				iteration += 8;
			}
			gatherLightPaths(m_pathSampler, m_config.useVC, m_config.useVM, radius, hilbertCurve.getPointCount(), result, batres);

			for (size_t i = 0; i < hilbertCurve.getPointCount(); ++i) {
				if (stop) break;

				Point2i offset = Point2i(hilbertCurve[i]);
				m_sampler->generate(offset);
				sampleCameraPath(m_pathSampler, result, m_config.useVC, m_config.useVM, radius, offset, i, *splats);

				for (size_t k = 0; k < splats->size(); ++k) {
					Spectrum value = splats->getValue(k);
					result->putSample(splats->getPosition(k), &value[0]);
					// [UC] for unbiased check					
					if (batres != NULL)
						batres->put(splats->getPosition(k), &value[0]);
				}
			}
			actualSampleCount++;

#if UPM_DEBUG == 1
			// [UC] for unbiased check
			if (m_config.enableSeparateDump){
				numSepSamples++;
				if (intervalTimer->getSeconds() >= sepInterval){
					Bitmap *bitmap = const_cast<Bitmap *>(batres->getBitmap());
					ref<Bitmap> hdrBitmap = bitmap->convert(Bitmap::ERGB, Bitmap::EFloat32, 1.0, 1.f / (Float)numSepSamples);
					fs::path filename = fs::path(formatString("E://%s_k%02d.pfm", "vcm_sep", numBatch * 8 + workID));
					ref<FileStream> targetFile = new FileStream(filename,
						FileStream::ETruncReadWrite);
					hdrBitmap->write(Bitmap::EPFM, targetFile, 1);
					batres->clear();
					intervalTimer->reset();
					numSepSamples = 0;
					numBatch++;
				}
			}

			// [UC] for relative contribution graphing
			if (m_config.enableProgressiveDump && intervalTimer2->getSeconds() >= rcInterval){
				fs::path path = m_scene->getDestinationFile();
				result->progressiveDump(filmSize.x, filmSize.y, m_config.maxDepth, path.parent_path(), path.stem(), numProgressiveBatch * 8 + workID, actualSampleCount, false);
				intervalTimer2->reset();
				numProgressiveBatch++;
			}
#endif
		}

#if UPM_DEBUG == 1
		// [UC] for relative contribution graphing
		if (m_config.enableProgressiveDump){
			fs::path path = m_scene->getDestinationFile();
			result->progressiveDump(filmSize.x, filmSize.y, m_config.maxDepth, path.parent_path(), path.stem(), numProgressiveBatch * 8 + workID, actualSampleCount, false);
			intervalTimer2->reset();
			numProgressiveBatch++;
		}
#endif

		Log(EInfo, "Run %d iterations", actualSampleCount);
		result->accumSampleCount(actualSampleCount);

		delete splats;
	}

	ref<WorkProcessor> clone() const {
		return new VCMRenderer(m_config);
	}

	MTS_DECLARE_CLASS()
private:
	VCMConfiguration m_config;
	ref<Scene> m_scene;
	ref<Sensor> m_sensor;
	ref<Film> m_film;
	ref<PathSampler> m_pathSampler;
	ref<Sampler> m_sampler;

	size_t m_lightPathNum;
	LightPathTreeV m_lightPathTree;
	std::vector<LightVertexV> m_lightVertices;
	std::vector<LightVertexExtV> m_lightVerticesExt;
	std::vector<size_t> m_lightPathEnds;
};

/* ==================================================================== */
/*                           Parallel process                           */
/* ==================================================================== */

VCMProcess::VCMProcess(const RenderJob *parent, RenderQueue *queue,
	const VCMConfiguration &conf) : m_job(parent), m_queue(queue),
		m_config(conf), m_progress(NULL) {
	m_timeoutTimer = new Timer();
	m_refreshTimer = new Timer();
	m_resultMutex = new Mutex();
	m_resultCounter = 0;
	m_workCounter = 0;
	m_refreshTimeout = 1;
}

ref<WorkProcessor> VCMProcess::createWorkProcessor() const {
	return new VCMRenderer(m_config);
}

void VCMProcess::develop() {
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

void VCMProcess::processResult(const WorkResult *workResult, bool cancelled) {
	const UPMWorkResult *wr = static_cast<const UPMWorkResult *>(workResult);
	LockGuard lock(m_resultMutex);	
	//m_accum->put(result);
	m_result->put(wr);
	m_progress->update(++m_resultCounter);
	m_refreshTimeout = std::min(2000U, m_refreshTimeout * 2);

	/* Re-develop the entire image every two seconds if partial results are
	   visible (e.g. in a graphical user interface). */
	if (m_job->isInteractive() && m_refreshTimer->getMilliseconds() > m_refreshTimeout)
		develop();
}

ParallelProcess::EStatus VCMProcess::generateWork(WorkUnit *unit, int worker) {
	int timeout = 0;
	if (m_config.timeout > 0) {
		timeout = static_cast<int>(static_cast<int64_t>(m_config.timeout*1000) -
		          static_cast<int64_t>(m_timeoutTimer->getMilliseconds()));
	}

	if (m_workCounter >= m_config.workUnits || timeout < 0)
		return EFailure;

	SeedWorkUnit *workUnit = static_cast<SeedWorkUnit *>(unit);
	workUnit->setID(m_workCounter++);
	workUnit->setTimeout(timeout);
	return ESuccess;
}

void VCMProcess::bindResource(const std::string &name, int id) {
	ParallelProcess::bindResource(name, id);
	if (name == "sensor") {
		m_film = static_cast<Sensor *>(Scheduler::getInstance()->getResource(id))->getFilm();
		if (m_progress)
			delete m_progress;
		m_progress = new ProgressReporter("Rendering", m_config.workUnits, m_job);
// 		m_accum = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize());
// 		m_accum->clear();
		m_result = new UPMWorkResult(m_film->getCropSize().x, m_film->getCropSize().y, m_config.maxDepth, NULL);
		m_result->clear();
		m_developBuffer = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, m_film->getCropSize());
	}
}

MTS_IMPLEMENT_CLASS_S(VCMRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(VCMProcess, false, ParallelProcess)
MTS_IMPLEMENT_CLASS(SeedWorkUnit, false, WorkUnit)
MTS_IMPLEMENT_CLASS(UPMWorkResult, false, WorkResult)
MTS_NAMESPACE_END
