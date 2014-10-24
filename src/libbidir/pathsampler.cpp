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

#include <mitsuba/bidir/pathsampler.h>
#include <mitsuba/bidir/rsampler.h>
#include <mitsuba/bidir/util.h>
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include <boost/bind.hpp>

MTS_NAMESPACE_BEGIN

PathSampler::PathSampler(ETechnique technique, const Scene *scene, Sampler *sensorSampler,
		Sampler *emitterSampler, Sampler *directSampler, int maxDepth, int rrDepth,
		bool excludeDirectIllum,  bool sampleDirect, bool lightImage,
		Float gatherRadius, Sampler *lightPathSampler)
	: m_technique(technique), m_scene(scene), m_emitterSampler(emitterSampler),
	  m_sensorSampler(sensorSampler), m_directSampler(directSampler), m_maxDepth(maxDepth),
	  m_rrDepth(rrDepth), m_excludeDirectIllum(excludeDirectIllum), m_sampleDirect(sampleDirect),
      m_lightImage(lightImage),
	  m_lightPathSampler(lightPathSampler){

	if (technique == EUnidirectional) {
		/* Instantiate a volumetric path tracer */
		Properties props("volpath");
		props.setInteger("maxDepth", maxDepth);
		props.setInteger("rrDepth", rrDepth);
		m_integrator = static_cast<SamplingIntegrator *> (PluginManager::getInstance()->
			createObject(MTS_CLASS(SamplingIntegrator), props));
	}

	/* Determine the necessary random walk depths based on properties of the endpoints */
	m_emitterDepth = m_sensorDepth = maxDepth;

	/* Go one extra step if the sensor can be intersected */
	if (!m_scene->hasDegenerateSensor() && m_emitterDepth != -1)
		++m_emitterDepth;

	/* Go one extra step if there are emitters that can be intersected */
	if (!m_scene->hasDegenerateEmitters() && m_sensorDepth != -1)
		++m_sensorDepth;
}

PathSampler::~PathSampler() {
	if (!m_pool.unused())
		Log(EWarn, "Warning: memory pool still contains used objects!");
}

void PathSampler::sampleSplats(const Point2i &offset, SplatList &list) {
	const Sensor *sensor = m_scene->getSensor();

	list.clear();
	switch (m_technique) {
		case EBidirectional: {
				/* Uniformly sample a scene time */
				Float time = sensor->getShutterOpen();
				if (sensor->needsTimeSample())
					time = sensor->sampleTime(m_sensorSampler->next1D());

				/* Initialize the path endpoints */
				m_emitterSubpath.initialize(m_scene, time, EImportance, m_pool);
				m_sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);

				/* Perform two random walks from the sensor and emitter side */
				m_emitterSubpath.randomWalk(m_scene, m_emitterSampler, m_emitterDepth,
						m_rrDepth, EImportance, m_pool);

				if (offset == Point2i(-1))
					m_sensorSubpath.randomWalk(m_scene, m_sensorSampler,
						m_sensorDepth, m_rrDepth, ERadiance, m_pool);
				else
					m_sensorSubpath.randomWalkFromPixel(m_scene, m_sensorSampler,
						m_sensorDepth, offset, m_rrDepth, m_pool);

				/* Compute the combined weights along the two subpaths */
				Spectrum *importanceWeights = (Spectrum *) alloca(m_emitterSubpath.vertexCount() * sizeof(Spectrum)),
						 *radianceWeights  = (Spectrum *) alloca(m_sensorSubpath.vertexCount()  * sizeof(Spectrum));

				importanceWeights[0] = radianceWeights[0] = Spectrum(1.0f);
				for (size_t i=1; i<m_emitterSubpath.vertexCount(); ++i)
					importanceWeights[i] = importanceWeights[i-1] *
						m_emitterSubpath.vertex(i-1)->weight[EImportance] *
						m_emitterSubpath.vertex(i-1)->rrWeight *
						m_emitterSubpath.edge(i-1)->weight[EImportance];

				for (size_t i=1; i<m_sensorSubpath.vertexCount(); ++i)
					radianceWeights[i] = radianceWeights[i-1] *
						m_sensorSubpath.vertex(i-1)->weight[ERadiance] *
						m_sensorSubpath.vertex(i-1)->rrWeight *
						m_sensorSubpath.edge(i-1)->weight[ERadiance];

				if (m_sensorSubpath.vertexCount() > 2) {
					Point2 samplePos(0.0f);
					m_sensorSubpath.vertex(1)->getSamplePosition(m_sensorSubpath.vertex(2), samplePos);
					list.append(samplePos, Spectrum(0.0f));
				}

				PathVertex tempEndpoint, tempSample;
				PathEdge tempEdge, connectionEdge;
				Point2 samplePos(0.0f);

				for (int s = (int) m_emitterSubpath.vertexCount()-1; s >= 0; --s) {
					/* Determine the range of sensor vertices to be traversed,
					   while respecting the specified maximum path length */
					int minT = std::max(2-s, m_lightImage ? 0 : 2),
					    maxT = (int) m_sensorSubpath.vertexCount() - 1;
					if (m_maxDepth != -1)
						maxT = std::min(maxT, m_maxDepth + 1 - s);

					for (int t = maxT; t >= minT; --t) {
						PathVertex
							*vsPred = m_emitterSubpath.vertexOrNull(s-1),
							*vtPred = m_sensorSubpath.vertexOrNull(t-1),
							*vs = m_emitterSubpath.vertex(s),
							*vt = m_sensorSubpath.vertex(t);
						PathEdge
							*vsEdge = m_emitterSubpath.edgeOrNull(s-1),
							*vtEdge = m_sensorSubpath.edgeOrNull(t-1);

						RestoreMeasureHelper rmh0(vs), rmh1(vt);

						/* Will be set to true if direct sampling was used */
						bool sampleDirect = false;

						/* Number of edges of the combined subpaths */
						int depth = s + t - 1;

						/* Allowed remaining number of ENull vertices that can
						   be bridged via pathConnect (negative=arbitrarily many) */
						int remaining = m_maxDepth - depth;

						/* Will receive the path weight of the (s, t)-connection */
						Spectrum value;

						/* Account for the terms of the measurement contribution
						   function that are coupled to the connection endpoints */
						if (vs->isEmitterSupernode()) {
							/* If possible, convert 'vt' into an emitter sample */
							if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
								continue;

							value = radianceWeights[t] *
								vs->eval(m_scene, vsPred, vt, EImportance) *
								vt->eval(m_scene, vtPred, vs, ERadiance);
						} else if (vt->isSensorSupernode()) {
							/* If possible, convert 'vs' into an sensor sample */
							if (!vs->cast(m_scene, PathVertex::ESensorSample) || vs->isDegenerate())
								continue;

							/* Make note of the changed pixel sample position */
							if (!vs->getSamplePosition(vsPred, samplePos))
								continue;

							value = importanceWeights[s] *
								vs->eval(m_scene, vsPred, vt, EImportance) *
								vt->eval(m_scene, vtPred, vs, ERadiance);
						} else if (m_sampleDirect && ((t == 1 && s > 1) || (s == 1 && t > 1))) {
							/* s==1/t==1 path: use a direct sampling strategy if requested */
							if (s == 1) {
								if (vt->isDegenerate())
									continue;

								/* Generate a position on an emitter using direct sampling */
								value = radianceWeights[t] * vt->sampleDirect(m_scene, m_directSampler,
									&tempEndpoint, &tempEdge, &tempSample, EImportance);

								if (value.isZero())
									continue;
								vs = &tempSample; vsPred = &tempEndpoint; vsEdge = &tempEdge;
								value *= vt->eval(m_scene, vtPred, vs, ERadiance);
								vt->measure = EArea;
							} else {
								if (vs->isDegenerate())
									continue;
								/* Generate a position on the sensor using direct sampling */
								value = importanceWeights[s] * vs->sampleDirect(m_scene, m_directSampler,
									&tempEndpoint, &tempEdge, &tempSample, ERadiance);
								if (value.isZero())
									continue;
								vt = &tempSample; vtPred = &tempEndpoint; vtEdge = &tempEdge;
								value *= vs->eval(m_scene, vsPred, vt, EImportance);
								vs->measure = EArea;
							}

							sampleDirect = true;
						} else {
							/* Can't connect degenerate endpoints */
							if (vs->isDegenerate() || vt->isDegenerate())
								continue;

							value = importanceWeights[s] * radianceWeights[t] *
								vs->eval(m_scene, vsPred, vt, EImportance) *
								vt->eval(m_scene, vtPred, vs, ERadiance);

							/* Temporarily force vertex measure to EArea. Needed to
							   handle BSDFs with diffuse + specular components */
							vs->measure = vt->measure = EArea;
						}

						/* Attempt to connect the two endpoints, which could result in
						   the creation of additional vertices (index-matched boundaries etc.) */
						int interactions = remaining;
						if (value.isZero() || !connectionEdge.pathConnectAndCollapse(
								m_scene, vsEdge, vs, vt, vtEdge, interactions))
							continue;

						depth += interactions;

						if (m_excludeDirectIllum && depth <= 2)
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
								m_sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
							else
								m_emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
						}

						/* Compute the multiple importance sampling weight */
						value *= Path::miWeight(false, m_scene, m_emitterSubpath, &connectionEdge,
							m_sensorSubpath, s, t, m_sampleDirect, m_lightImage);

						if (sampleDirect) {
							/* Now undo the previous change */
							if (t == 1)
								m_sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
							else
								m_emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
						}

						/* Determine the pixel sample position when necessary */
						if (vt->isSensorSample() && !vt->getSamplePosition(vs, samplePos))
							continue;


						if (t < 2) {
							list.append(samplePos, value);
						} else {
							BDAssert(m_sensorSubpath.vertexCount() > 2);
							list.accum(0, value);
						}
					}
				}

				/* Release any used edges and vertices back to the memory pool */
				m_sensorSubpath.release(m_pool);
				m_emitterSubpath.release(m_pool);
			}
			break;

		case EUnidirectional: {
				Point2 apertureSample(0.5f);
				Float timeSample = 0.5f;

				if (sensor->needsApertureSample())
					apertureSample = m_sensorSampler->next2D();
				if (sensor->needsTimeSample())
					timeSample = m_sensorSampler->next1D();

				const Vector2i cropSize = sensor->getFilm()->getCropSize();
				Point2 samplePos = m_sensorSampler->next2D();
				if (offset == Point2i(-1)) {
					samplePos.x *= cropSize.x;
					samplePos.y *= cropSize.y;
				} else {
					samplePos.x += offset.x;
					samplePos.y += offset.y;
				}

				RayDifferential sensorRay;
				Spectrum value = sensor->sampleRayDifferential(
					sensorRay, samplePos, apertureSample, timeSample);

				RadianceQueryRecord rRec(m_scene, m_sensorSampler);
				rRec.newQuery(
					m_excludeDirectIllum ? (RadianceQueryRecord::ERadiance
					& ~(RadianceQueryRecord::EDirectSurfaceRadiance | RadianceQueryRecord::EEmittedRadiance))
					: RadianceQueryRecord::ERadiance, sensor->getMedium());

				value *= m_integrator->Li(sensorRay, rRec);

				list.append(samplePos, value);
			}
			break;
		default:
			Log(EError, "PathSampler::sample(): invalid technique!");
	}
}

void PathSampler::samplePaths(const Point2i &offset, PathCallback &callback) {
	BDAssert(m_technique == EBidirectional);

	const Sensor *sensor = m_scene->getSensor();
	/* Uniformly sample a scene time */
	Float time = sensor->getShutterOpen();
	if (sensor->needsTimeSample())
		time = sensor->sampleTime(m_sensorSampler->next1D());

	/* Initialize the path endpoints */
	m_emitterSubpath.initialize(m_scene, time, EImportance, m_pool);
	m_sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);

	/* Perform two random walks from the sensor and emitter side */
	m_emitterSubpath.randomWalk(m_scene, m_emitterSampler, m_emitterDepth,
			m_rrDepth, EImportance, m_pool);

	if (offset == Point2i(-1))
		m_sensorSubpath.randomWalk(m_scene, m_sensorSampler,
			m_sensorDepth, m_rrDepth, ERadiance, m_pool);
	else
		m_sensorSubpath.randomWalkFromPixel(m_scene, m_sensorSampler,
			m_sensorDepth, offset, m_rrDepth, m_pool);

	/* Compute the combined weights along the two subpaths */
	Spectrum *importanceWeights = (Spectrum *) alloca(m_emitterSubpath.vertexCount() * sizeof(Spectrum)),
			 *radianceWeights  = (Spectrum *) alloca(m_sensorSubpath.vertexCount()  * sizeof(Spectrum));

	importanceWeights[0] = radianceWeights[0] = Spectrum(1.0f);
	for (size_t i=1; i<m_emitterSubpath.vertexCount(); ++i)
		importanceWeights[i] = importanceWeights[i-1] *
			m_emitterSubpath.vertex(i-1)->weight[EImportance] *
			m_emitterSubpath.vertex(i-1)->rrWeight *
			m_emitterSubpath.edge(i-1)->weight[EImportance];

	for (size_t i=1; i<m_sensorSubpath.vertexCount(); ++i)
		radianceWeights[i] = radianceWeights[i-1] *
			m_sensorSubpath.vertex(i-1)->weight[ERadiance] *
			m_sensorSubpath.vertex(i-1)->rrWeight *
			m_sensorSubpath.edge(i-1)->weight[ERadiance];

	PathVertex tempEndpoint, tempSample, vsTemp, vtTemp;
	PathEdge tempEdge, connectionEdge;
	Point2 samplePos(0.0f);

	for (int s = (int) m_emitterSubpath.vertexCount()-1; s >= 0; --s) {
		/* Determine the range of sensor vertices to be traversed,
		   while respecting the specified maximum path length */
		int minT = std::max(2-s, m_lightImage ? 0 : 2),
		    maxT = (int) m_sensorSubpath.vertexCount() - 1;
		if (m_maxDepth != -1)
			maxT = std::min(maxT, m_maxDepth + 1 - s);

		for (int t = maxT; t >= minT; --t) {
			PathVertex
				*vsPred = m_emitterSubpath.vertexOrNull(s-1),
				*vtPred = m_sensorSubpath.vertexOrNull(t-1),
				*vs = m_emitterSubpath.vertex(s),
				*vt = m_sensorSubpath.vertex(t);
			PathEdge
				*vsEdge = m_emitterSubpath.edgeOrNull(s-1),
				*vtEdge = m_sensorSubpath.edgeOrNull(t-1);

			RestoreMeasureHelper rmh0(vs), rmh1(vt);

			/* Will be set to true if direct sampling was used */
			bool sampleDirect = false;

			/* Number of edges of the combined subpaths */
			int depth = s + t - 1;

			/* Allowed remaining number of ENull vertices that can
			   be bridged via pathConnect (negative=arbitrarily many) */
			int remaining = m_maxDepth - depth;

			/* Will receive the path weight of the (s, t)-connection */
			Spectrum value;

			/* Measure associated with the connection vertices */
			EMeasure vsMeasure = EArea, vtMeasure = EArea;

			/* Account for the terms of the measurement contribution
			   function that are coupled to the connection endpoints */
			if (vs->isEmitterSupernode()) {
				/* If possible, convert 'vt' into an emitter sample */
				if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
					continue;

				value = radianceWeights[t] *
					vs->eval(m_scene, vsPred, vt, EImportance) *
					vt->eval(m_scene, vtPred, vs, ERadiance);
			} else if (vt->isSensorSupernode()) {
				/* If possible, convert 'vs' into an sensor sample */
				if (!vs->cast(m_scene, PathVertex::ESensorSample) || vs->isDegenerate())
					continue;

				/* Make note of the changed pixel sample position */
				if (!vs->getSamplePosition(vsPred, samplePos))
					continue;

				value = importanceWeights[s] *
					vs->eval(m_scene, vsPred, vt, EImportance) *
					vt->eval(m_scene, vtPred, vs, ERadiance);
			} else if (m_sampleDirect && ((t == 1 && s > 1) || (s == 1 && t > 1))) {
				/* s==1/t==1 path: use a direct sampling strategy if requested */
				if (s == 1) {
					if (vt->isDegenerate())
						continue;

					/* Generate a position on an emitter using direct sampling */
					value = radianceWeights[t] * vt->sampleDirect(m_scene, m_directSampler,
						&tempEndpoint, &tempEdge, &tempSample, EImportance);

					if (value.isZero())
						continue;
					vs = &tempSample; vsPred = &tempEndpoint; vsEdge = &tempEdge;
					value *= vt->eval(m_scene, vtPred, vs, ERadiance);
					vsMeasure = vs->getAbstractEmitter()->needsDirectionSample() ? EArea : EDiscrete;
					vt->measure = EArea;
				} else {
					if (vs->isDegenerate())
						continue;
					/* Generate a position on the sensor using direct sampling */
					value = importanceWeights[s] * vs->sampleDirect(m_scene, m_directSampler,
						&tempEndpoint, &tempEdge, &tempSample, ERadiance);
					if (value.isZero())
						continue;
					vt = &tempSample; vtPred = &tempEndpoint; vtEdge = &tempEdge;
					value *= vs->eval(m_scene, vsPred, vt, EImportance);
					vtMeasure = vt->getAbstractEmitter()->needsDirectionSample() ? EArea : EDiscrete;
					vs->measure = EArea;
				}

				sampleDirect = true;
			} else {
				/* Can't connect degenerate endpoints */
				if (vs->isDegenerate() || vt->isDegenerate())
					continue;

				value = importanceWeights[s] * radianceWeights[t] *
					vs->eval(m_scene, vsPred, vt, EImportance) *
					vt->eval(m_scene, vtPred, vs, ERadiance);

				/* Temporarily force vertex measure to EArea. Needed to
				   handle BSDFs with diffuse + specular components */
				vs->measure = vt->measure = EArea;
			}

			/* Attempt to connect the two endpoints, which could result in
			   the creation of additional vertices (index-matched boundaries etc.) */
			m_connectionSubpath.release(m_pool);
			if (value.isZero() || !PathEdge::pathConnect(m_scene, vsEdge,
					vs, m_connectionSubpath, vt, vtEdge, remaining, m_pool))
				continue;

			depth += (int) m_connectionSubpath.vertexCount();

			if (m_excludeDirectIllum && depth <= 2)
				continue;

			PathEdge connectionEdge;
			m_connectionSubpath.collapseTo(connectionEdge);

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
					m_sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
				else
					m_emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
			}

			/* Compute the multiple importance sampling weight */
			value *= Path::miWeight(false, m_scene, m_emitterSubpath, &connectionEdge,
				m_sensorSubpath, s, t, m_sampleDirect, m_lightImage);

			if (!value.isZero()) {
				int k = (int) m_connectionSubpath.vertexCount();
				/* Construct the full path, make a temporary backup copy of the connection vertices */
				m_fullPath.clear();
				m_fullPath.append(m_emitterSubpath, 0, s+1);
				m_fullPath.append(m_connectionSubpath);
				m_fullPath.append(m_sensorSubpath, 0, t+1, true);
				m_fullPath.vertex(s)     = &vsTemp;
				m_fullPath.vertex(s+k+1) = &vtTemp;
				vsTemp = *m_emitterSubpath.vertex(s);
				vtTemp = *m_sensorSubpath.vertex(t);

				if (vsTemp.update(m_scene, m_fullPath.vertexOrNull(s-1),   m_fullPath.vertex(s+1), EImportance, vsMeasure) &&
				    vtTemp.update(m_scene, m_fullPath.vertexOrNull(s+k+2), m_fullPath.vertex(s+k), ERadiance, vtMeasure)) {

					if (vtTemp.isSensorSample())
						vtTemp.updateSamplePosition(&vsTemp);

					callback(s, t, value.getLuminance(), m_fullPath);
				} else {
					Log(EWarn, "samplePaths(): internal error in update() "
						"failed! : %s, %s", vs->toString().c_str(), vt->toString().c_str());
				}
			}

			if (sampleDirect) {
				/* Now undo the previous change */
				if (t == 1)
					m_sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
				else
					m_emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
			}
		}
	}

	/* Release any used edges and vertices back to the memory pool */
	m_sensorSubpath.release(m_pool);
	m_emitterSubpath.release(m_pool);
	m_connectionSubpath.release(m_pool);
}

/**
 * \brief Sort predicate used to order \ref PathSeed instances by
 * their index into the underlying random number stream
 */
struct PathSeedSortPredicate {
	bool operator()(const PathSeed &left, const PathSeed &right) {
		return left.sampleIndex < right.sampleIndex;
	}
};

Float PathSampler::computeAverageLuminance(size_t sampleCount) {
	Log(EInfo, "Integrating luminance values over the image plane ("
			SIZE_T_FMT " samples)..", sampleCount);

	ref<Timer> timer = new Timer();

	SplatList splatList;
	Float mean = 0.0f, variance = 0.0f;
	for (size_t i=0; i<sampleCount; ++i) {
		/* Run the path sampling strategy */
		m_sensorSampler->generate(Point2i(0));
		sampleSplats(Point2i(-1), splatList);

		Float lum = splatList.luminance,
		      delta = lum - mean;
		mean += delta / (Float) (i+1);
		variance += delta * (lum - mean);
	}

	BDAssert(m_pool.unused());
	Float stddev = std::sqrt(variance / (sampleCount-1));

	Log(EInfo, "Done -- average luminance value = %f, stddev = %f (took %i ms)",
			mean, stddev, timer->getMilliseconds());

	if (mean == 0)
		Log(EError, "The average image luminance appears to be zero! This could indicate "
			"a problem with the scene setup. Aborting the rendering process.");

	return mean;
}

static void seedCallback(std::vector<PathSeed> &output, const Bitmap *importanceMap,
		Float &accum, int s, int t, Float weight, Path &path) {
	accum += weight;

	if (importanceMap) {
		const Float *luminanceValues = importanceMap->getFloatData();
		Vector2i size = importanceMap->getSize();

		const Point2 &pos = path.getSamplePosition();
		Point2i intPos(
			std::min(std::max(0, (int) pos.x), size.x-1),
			std::min(std::max(0, (int) pos.y), size.y-1));
		weight /= luminanceValues[intPos.x + intPos.y * size.x];
	}

	output.push_back(PathSeed(0, weight, s, t));
}

Float PathSampler::generateSeeds(size_t sampleCount, size_t seedCount,
		bool fineGrained, const Bitmap *importanceMap, std::vector<PathSeed> &seeds) {
	Log(EInfo, "Integrating luminance values over the image plane ("
			SIZE_T_FMT " samples)..", sampleCount);

	BDAssert(m_sensorSampler == m_emitterSampler);
	BDAssert(m_sensorSampler->getClass()->derivesFrom(MTS_CLASS(ReplayableSampler)));

	ref<Timer> timer = new Timer();
	std::vector<PathSeed> tempSeeds;
	tempSeeds.reserve(sampleCount);

	SplatList splatList;
	Float luminance;
	PathCallback callback = boost::bind(&seedCallback,
		boost::ref(tempSeeds), importanceMap, boost::ref(luminance),
		_1, _2, _3, _4);

	Float mean = 0.0f, variance = 0.0f;
	for (size_t i=0; i<sampleCount; ++i) {
		size_t seedIndex = tempSeeds.size();
		size_t sampleIndex = m_sensorSampler->getSampleIndex();
		luminance = 0.0f;

		if (fineGrained) {
			samplePaths(Point2i(-1), callback);

			/* Fine seed granularity (e.g. for Veach-MLT).
			   Set the correct the sample index value */
			for (size_t j = seedIndex; j<tempSeeds.size(); ++j)
				tempSeeds[j].sampleIndex = sampleIndex;
		} else {
			/* Run the path sampling strategy */
			sampleSplats(Point2i(-1), splatList);
			luminance = splatList.luminance;
			splatList.normalize(importanceMap);

			/* Coarse seed granularity (e.g. for PSSMLT) */
			if (luminance != 0)
				tempSeeds.push_back(PathSeed(sampleIndex, luminance));
		}

		/* Numerically robust online variance estimation using an
		   algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232) */
		Float delta = luminance - mean;
		mean += delta / (Float) (i+1);
		variance += delta * (luminance - mean);
	}
	BDAssert(m_pool.unused());
	Float stddev = std::sqrt(variance / (sampleCount-1));

	Log(EInfo, "Done -- average luminance value = %f, stddev = %f (took %i ms)",
			mean, stddev, timer->getMilliseconds());

	if (mean == 0)
		Log(EError, "The average image luminance appears to be zero! This could indicate "
			"a problem with the scene setup. Aborting the MLT rendering process.");

	Log(EDebug, "Sampling " SIZE_T_FMT "/" SIZE_T_FMT " MLT seeds",
		seedCount, tempSeeds.size());

	DiscreteDistribution seedPDF(tempSeeds.size());
	for (size_t i=0; i<tempSeeds.size(); ++i)
		seedPDF.append(tempSeeds[i].luminance);
	seedPDF.normalize();

	seeds.clear();
	seeds.reserve(seedCount);
	for (size_t i=0; i<seedCount; ++i)
		seeds.push_back(tempSeeds.at(seedPDF.sample(m_sensorSampler->next1D())));

	/* Sort the seeds to avoid unnecessary rewinds in the ReplayableSampler */
	std::sort(seeds.begin(), seeds.end(), PathSeedSortPredicate());

	return mean;
}

static void reconstructCallback(const PathSeed &seed, const Bitmap *importanceMap,
		Path &result, MemoryPool &pool, int s, int t, Float weight, Path &path) {
	if (s == seed.s && t == seed.t) {
		if (importanceMap) {
			const Float *luminanceValues = importanceMap->getFloatData();
			Vector2i size = importanceMap->getSize();

			const Point2 &pos = path.getSamplePosition();
			Point2i intPos(
				std::min(std::max(0, (int) pos.x), size.x-1),
				std::min(std::max(0, (int) pos.y), size.y-1));
			weight /= luminanceValues[intPos.x + intPos.y * size.x];
		}

		if (seed.luminance != weight)
			SLog(EError, "Internal error in reconstructPath(): luminances "
				"don't match (%f vs %f)!", weight, seed.luminance);
		path.clone(result, pool);
	}
}

void PathSampler::reconstructPath(const PathSeed &seed, const Bitmap *importanceMap, Path &result) {
	ReplayableSampler *rplSampler = static_cast<ReplayableSampler *>(m_sensorSampler.get());
	Assert(result.length() == 0);

	/* Generate the initial sample by replaying the seeding random
	   number stream at the appropriate position. */
	rplSampler->setSampleIndex(seed.sampleIndex);

	PathCallback callback = boost::bind(&reconstructCallback,
		boost::cref(seed), importanceMap,
		boost::ref(result), boost::ref(m_pool), _1, _2, _3, _4);

	samplePaths(Point2i(-1), callback);

	if (result.length() == 0)
		Log(EError, "Internal error in reconstructPath(): desired configuration was never created!");
}

void SplatList::normalize(const Bitmap *importanceMap) {
	if (importanceMap) {
		luminance = 0.0f;

		/* Two-stage MLT -- weight contributions using a luminance image */
		const Float *luminanceValues = importanceMap->getFloatData();
		Vector2i size = importanceMap->getSize();
		for (size_t i=0; i<splats.size(); ++i) {
			if (splats[i].second.isZero())
				continue;

			const Point2 &pos = splats[i].first;
			Point2i intPos(
				std::min(std::max(0, (int) pos.x), size.x-1),
				std::min(std::max(0, (int) pos.y), size.y-1));
			Float lumValue = luminanceValues[intPos.x + intPos.y * size.x];
			splats[i].second /= lumValue;
			luminance += splats[i].second.getLuminance();
		}
	}

	if (luminance > 0) {
		/* Normalize the contributions */
		Float invLuminance = 1.0f / luminance;
		for (size_t i=0; i<splats.size(); ++i)
			splats[i].second *= invLuminance;
	}
}

/* Specific copy of functions for Connection MLT*/
void PathSampler::sampleSplatsConnection(const Point2i &offset, SplatListImp &list, const int connectionFlag) {
	const Sensor *sensor = m_scene->getSensor();

	list.clear();
	switch (m_technique) {
	case EBidirectional: {
		/* Uniformly sample a scene time */
		Float time = sensor->getShutterOpen();
		if (sensor->needsTimeSample())
			time = sensor->sampleTime(m_sensorSampler->next1D());

		/* Initialize the path endpoints */
		m_emitterSubpath.initialize(m_scene, time, EImportance, m_pool);
		m_sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);

		/* Perform two random walks from the sensor and emitter side */
		m_emitterSubpath.randomWalk(m_scene, m_emitterSampler, m_emitterDepth,
			m_rrDepth, EImportance, m_pool);

		if (offset == Point2i(-1))
			m_sensorSubpath.randomWalk(m_scene, m_sensorSampler,
			m_sensorDepth, m_rrDepth, ERadiance, m_pool);
		else
			m_sensorSubpath.randomWalkFromPixel(m_scene, m_sensorSampler,
			m_sensorDepth, offset, m_rrDepth, m_pool);

		/* Compute the combined weights along the two subpaths */
		Spectrum *importanceWeights = (Spectrum *)alloca(m_emitterSubpath.vertexCount() * sizeof(Spectrum)),
			*radianceWeights = (Spectrum *)alloca(m_sensorSubpath.vertexCount()  * sizeof(Spectrum));

		importanceWeights[0] = radianceWeights[0] = Spectrum(1.0f);
		for (size_t i = 1; i<m_emitterSubpath.vertexCount(); ++i)
			importanceWeights[i] = importanceWeights[i - 1] *
			m_emitterSubpath.vertex(i - 1)->weight[EImportance] *
			m_emitterSubpath.vertex(i - 1)->rrWeight *
			m_emitterSubpath.edge(i - 1)->weight[EImportance];

		for (size_t i = 1; i<m_sensorSubpath.vertexCount(); ++i)
			radianceWeights[i] = radianceWeights[i - 1] *
			m_sensorSubpath.vertex(i - 1)->weight[ERadiance] *
			m_sensorSubpath.vertex(i - 1)->rrWeight *
			m_sensorSubpath.edge(i - 1)->weight[ERadiance];

		if (m_sensorSubpath.vertexCount() > 2) {
			Point2 samplePos(0.0f);
			m_sensorSubpath.vertex(1)->getSamplePosition(m_sensorSubpath.vertex(2), samplePos);
			list.append(samplePos, Spectrum(0.0f), Spectrum(0.0f));
		}

		PathVertex tempEndpoint, tempSample;
		PathEdge tempEdge, connectionEdge;
		Point2 samplePos(0.0f);

		for (int s = (int)m_emitterSubpath.vertexCount() - 1; s >= 0; --s) {
			/* Determine the range of sensor vertices to be traversed,
			while respecting the specified maximum path length */
			int minT = std::max(2 - s, m_lightImage ? 0 : 2),
				maxT = (int)m_sensorSubpath.vertexCount() - 1;
			if (m_maxDepth != -1)
				maxT = std::min(maxT, m_maxDepth + 1 - s);

			for (int t = maxT; t >= minT; --t) {
				PathVertex
					*vsPred = m_emitterSubpath.vertexOrNull(s - 1),
					*vtPred = m_sensorSubpath.vertexOrNull(t - 1),
					*vs = m_emitterSubpath.vertex(s),
					*vt = m_sensorSubpath.vertex(t);
				PathEdge
					*vsEdge = m_emitterSubpath.edgeOrNull(s - 1),
					*vtEdge = m_sensorSubpath.edgeOrNull(t - 1);

				RestoreMeasureHelper rmh0(vs), rmh1(vt);

				/* Will be set to true if direct sampling was used */
				bool sampleDirect = false;

				/* Number of edges of the combined subpaths */
				int depth = s + t - 1;

				/* Allowed remaining number of ENull vertices that can
				be bridged via pathConnect (negative=arbitrarily many) */
				int remaining = m_maxDepth - depth;

				/* Will receive the path weight of the (s, t)-connection */
				Spectrum value;

				Spectrum importance = Spectrum(1.f);

				/* Account for the terms of the measurement contribution
				function that are coupled to the connection endpoints */
				if (vs->isEmitterSupernode()) {
					/* If possible, convert 'vt' into an emitter sample */
					if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
						continue;

					value = radianceWeights[t] *
						vs->eval(m_scene, vsPred, vt, EImportance) *
						vt->eval(m_scene, vtPred, vs, ERadiance);

					if (connectionFlag & EConnectImportance) importance *= radianceWeights[t];
					if (connectionFlag & EConnectBRDF) importance *= vs->eval(m_scene, vsPred, vt, EImportance) *	vt->eval(m_scene, vtPred, vs, ERadiance);					
				}
				else if (vt->isSensorSupernode()) {
					/* If possible, convert 'vs' into an sensor sample */
					if (!vs->cast(m_scene, PathVertex::ESensorSample) || vs->isDegenerate())
						continue;

					/* Make note of the changed pixel sample position */
					if (!vs->getSamplePosition(vsPred, samplePos))
						continue;

					value = importanceWeights[s] *
						vs->eval(m_scene, vsPred, vt, EImportance) *
						vt->eval(m_scene, vtPred, vs, ERadiance);

					if (connectionFlag & EConnectRadiance) importance *= importanceWeights[t];
					if (connectionFlag & EConnectBRDF) importance *= vs->eval(m_scene, vsPred, vt, EImportance) * 	vt->eval(m_scene, vtPred, vs, ERadiance);
				}
				else if (m_sampleDirect && ((t == 1 && s > 1) || (s == 1 && t > 1))) {
					/* s==1/t==1 path: use a direct sampling strategy if requested */
					if (s == 1) {
						if (vt->isDegenerate())
							continue;

						/* Generate a position on an emitter using direct sampling */
						value = radianceWeights[t] * vt->sampleDirect(m_scene, m_directSampler,
							&tempEndpoint, &tempEdge, &tempSample, EImportance);

						if (value.isZero())
							continue;
						vs = &tempSample; vsPred = &tempEndpoint; vsEdge = &tempEdge;
						value *= vt->eval(m_scene, vtPred, vs, ERadiance);
						vt->measure = EArea;

						if (connectionFlag & EConnectImportance) importance *= radianceWeights[t];
						if (connectionFlag & EConnectRadiance) importance *= vt->sampleDirect(m_scene, m_directSampler, &tempEndpoint, &tempEdge, &tempSample, EImportance);
						if (connectionFlag & EConnectBRDF) importance *= vt->eval(m_scene, vtPred, vs, ERadiance);
					}
					else {
						if (vs->isDegenerate())
							continue;
						/* Generate a position on the sensor using direct sampling */
						value = importanceWeights[s] * vs->sampleDirect(m_scene, m_directSampler,
							&tempEndpoint, &tempEdge, &tempSample, ERadiance);
						if (value.isZero())
							continue;
						vt = &tempSample; vtPred = &tempEndpoint; vtEdge = &tempEdge;
						value *= vs->eval(m_scene, vsPred, vt, EImportance);
						vs->measure = EArea;

						if (connectionFlag & EConnectImportance) importance *= vs->sampleDirect(m_scene, m_directSampler, &tempEndpoint, &tempEdge, &tempSample, ERadiance);
						if (connectionFlag & EConnectRadiance) importance *= importanceWeights[s];
						if (connectionFlag & EConnectBRDF) importance *= vs->eval(m_scene, vsPred, vt, EImportance);
					}

					sampleDirect = true;
				}
				else {
					/* Can't connect degenerate endpoints */
					if (vs->isDegenerate() || vt->isDegenerate())
						continue;

					value = importanceWeights[s] * radianceWeights[t] *
						vs->eval(m_scene, vsPred, vt, EImportance) *
						vt->eval(m_scene, vtPred, vs, ERadiance);

					if (connectionFlag & EConnectImportance) importance *= radianceWeights[t];
					if (connectionFlag & EConnectRadiance) importance *= importanceWeights[s];
					if (connectionFlag & EConnectBRDF) importance *= vs->eval(m_scene, vsPred, vt, EImportance) *	vt->eval(m_scene, vtPred, vs, ERadiance);

					/* Temporarily force vertex measure to EArea. Needed to
					handle BSDFs with diffuse + specular components */
					vs->measure = vt->measure = EArea;
				}

				/* Attempt to connect the two endpoints, which could result in
				the creation of additional vertices (index-matched boundaries etc.) */
				int interactions = remaining;
				if (value.isZero() || !connectionEdge.pathConnectAndCollapse(
					m_scene, vsEdge, vs, vt, vtEdge, interactions))
					continue;

				depth += interactions;

				if (m_excludeDirectIllum && depth <= 2)
					continue;

				/* Account for the terms of the measurement contribution
				function that are coupled to the connection edge */
				Spectrum interThroughput = Spectrum(0.0f);
				if (!sampleDirect)
					interThroughput = connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
				else
					interThroughput = connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
					(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));
				value *= interThroughput;				

				if (sampleDirect) {
					/* A direct sampling strategy was used, which generated
					two new vertices at one of the path ends. Temporarily
					modify the path to reflect this change */
					if (t == 1)
						m_sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
					else
						m_emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
				}

				/* Compute the multiple importance sampling weight */
				Float misWeight = Path::miWeight(false, m_scene, m_emitterSubpath, &connectionEdge,
					m_sensorSubpath, s, t, m_sampleDirect, m_lightImage);
				value *= misWeight;

				if (connectionFlag & EConnectGeometry)	importance *= interThroughput;
				if (connectionFlag & EConnectMis) importance *= misWeight;
				if (connectionFlag & EConnectAll) importance = value;


				if (sampleDirect) {
					/* Now undo the previous change */
					if (t == 1)
						m_sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
					else
						m_emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
				}

				/* Determine the pixel sample position when necessary */
				if (vt->isSensorSample() && !vt->getSamplePosition(vs, samplePos))
					continue;


				if (t < 2) {
					list.append(samplePos, value, importance);
				}
				else {
					BDAssert(m_sensorSubpath.vertexCount() > 2);
					list.accum(0, value, importance);
				}
			}
		}

		/* Release any used edges and vertices back to the memory pool */
		m_sensorSubpath.release(m_pool);
		m_emitterSubpath.release(m_pool);
	}
		break;

	case EUnidirectional: {
		Point2 apertureSample(0.5f);
		Float timeSample = 0.5f;

		if (sensor->needsApertureSample())
			apertureSample = m_sensorSampler->next2D();
		if (sensor->needsTimeSample())
			timeSample = m_sensorSampler->next1D();

		const Vector2i cropSize = sensor->getFilm()->getCropSize();
		Point2 samplePos = m_sensorSampler->next2D();
		if (offset == Point2i(-1)) {
			samplePos.x *= cropSize.x;
			samplePos.y *= cropSize.y;
		}
		else {
			samplePos.x += offset.x;
			samplePos.y += offset.y;
		}

		RayDifferential sensorRay;
		Spectrum value = sensor->sampleRayDifferential(
			sensorRay, samplePos, apertureSample, timeSample);

		RadianceQueryRecord rRec(m_scene, m_sensorSampler);
		rRec.newQuery(
			m_excludeDirectIllum ? (RadianceQueryRecord::ERadiance
			& ~(RadianceQueryRecord::EDirectSurfaceRadiance | RadianceQueryRecord::EEmittedRadiance))
			: RadianceQueryRecord::ERadiance, sensor->getMedium());

		value *= m_integrator->Li(sensorRay, rRec);

		list.append(samplePos, value, value);
	}
		break;
	default:
		Log(EError, "PathSampler::sample(): invalid technique!");
	}
}
Float PathSampler::generateSeedsConnection(size_t sampleCount, size_t seedCount,
	bool fineGrained, const Bitmap *importanceMap, std::vector<PathSeed> &seeds,
	const int connectionFlag) {
	Log(EInfo, "Integrating luminance values over the image plane ("
		SIZE_T_FMT " samples)..", sampleCount);

	BDAssert(m_sensorSampler == m_emitterSampler);
	BDAssert(m_sensorSampler->getClass()->derivesFrom(MTS_CLASS(ReplayableSampler)));

	ref<Timer> timer = new Timer();
	std::vector<PathSeed> tempSeeds;
	tempSeeds.reserve(sampleCount);

	SplatListImp splatList;
	Float luminance;
	PathCallback callback = boost::bind(&seedCallback,
		boost::ref(tempSeeds), importanceMap, boost::ref(luminance),
		_1, _2, _3, _4);

	Float mean = 0.0f, variance = 0.0f;
	for (size_t i = 0; i<sampleCount; ++i) {
		size_t seedIndex = tempSeeds.size();
		size_t sampleIndex = m_sensorSampler->getSampleIndex();
		luminance = 0.0f;

		if (fineGrained) {
			samplePaths(Point2i(-1), callback);

			/* Fine seed granularity (e.g. for Veach-MLT).
			Set the correct the sample index value */
			for (size_t j = seedIndex; j<tempSeeds.size(); ++j)
				tempSeeds[j].sampleIndex = sampleIndex;
		}
		else {
			/* Run the path sampling strategy */
			sampleSplatsConnection(Point2i(-1), splatList, connectionFlag);
			luminance = splatList.importance;
			splatList.normalize(importanceMap);

			/* Coarse seed granularity (e.g. for PSSMLT) */
			if (luminance != 0)
				tempSeeds.push_back(PathSeed(sampleIndex, luminance));
		}

		/* Numerically robust online variance estimation using an
		algorithm proposed by Donald Knuth (TAOCP vol.2, 3rd ed., p.232) */
		Float delta = luminance - mean;
		mean += delta / (Float)(i + 1);
		variance += delta * (luminance - mean);
	}
	BDAssert(m_pool.unused());
	Float stddev = std::sqrt(variance / (sampleCount - 1));

	Log(EInfo, "Done -- average luminance value = %f, stddev = %f (took %i ms)",
		mean, stddev, timer->getMilliseconds());

	if (mean == 0)
		Log(EError, "The average image luminance appears to be zero! This could indicate "
		"a problem with the scene setup. Aborting the MLT rendering process.");

	Log(EDebug, "Sampling " SIZE_T_FMT "/" SIZE_T_FMT " MLT seeds",
		seedCount, tempSeeds.size());

	DiscreteDistribution seedPDF(tempSeeds.size());
	for (size_t i = 0; i<tempSeeds.size(); ++i)
		seedPDF.append(tempSeeds[i].luminance);
	seedPDF.normalize();

	seeds.clear();
	seeds.reserve(seedCount);
	for (size_t i = 0; i<seedCount; ++i)
		seeds.push_back(tempSeeds.at(seedPDF.sample(m_sensorSampler->next1D())));

	/* Sort the seeds to avoid unnecessary rewinds in the ReplayableSampler */
	std::sort(seeds.begin(), seeds.end(), PathSeedSortPredicate());

	return mean;
}


/*
*	Tracing kernel for VCM
*/
void updateMisHelper(int i, const Path &path, MisState &state, const Scene* scene,
	size_t nLightPaths, Float misVcWeightFactor, Float misVmWeightFactor, ETransportMode mode){
	if (i == 0){
		state[EVCM] = 0.f;
		state[EVC] = 0.f;
		state[EVM] = 0.f;
	}
	else if (i == 1){
		state[EVCM] = 0.f;
		state[EVC] = 0.f;
		state[EVM] = 0.f;
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
			state[EVCM] = pconnect / (ptrace * p1);
			state[EVC] = 0.f;
			state[EVM] = 0.f;
		}
		else{
			state[EVCM] = pconnect / (ptrace * p1);
			if (measure != EDiscrete){
				Float geoTerm0 = v->isOnSurface() ? std::abs(dot(e->d, v->getGeometricNormal()) / (e->length * e->length)) : 1;
				state[EVC] = geoTerm0 / (ptrace * p1);
			}
			else{
				state[EVC] = 0.f;
			}
			state[EVM] = state[EVC] * misVcWeightFactor;
		}
	}
	else{
		PathVertex *v = path.vertex(i);
		PathEdge   *e = path.edge(i);
		if (v->measure == EDiscrete){
			PathVertex *vNext = path.vertex(i + 1);
			state[EVCM] = 0.f;
			state[EVC] *= std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / std::abs(vNext->isOnSurface() ? dot(e->d, vNext->getGeometricNormal()) : 1);
			state[EVM] *= std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / std::abs(vNext->isOnSurface() ? dot(e->d, vNext->getGeometricNormal()) : 1);
		}
		else{
			BDAssert(v->measure == EArea);
			PathVertex *vPred = path.vertex(i - 1);
			PathEdge * ePred = path.edge(i - 1);
			Float giIn = std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / (e->length * e->length);
			Float giInPred = std::abs(vPred->isOnSurface() ? dot(ePred->d, vPred->getGeometricNormal()) : 1) / (ePred->length * ePred->length);
			Float pi = v->pdf[mode] * e->pdf[mode];
			Float pir2_w = v->pdf[1 - mode] * ePred->pdf[1 - mode] / giInPred;
			state[EVC] = giIn / pi * (state[EVCM] + misVmWeightFactor + pir2_w * state[EVC]);
			state[EVM] = giIn / pi * (1.f + state[EVCM] * misVcWeightFactor + pir2_w * state[EVM]);
			state[EVCM] = 1 / pi;
		}
	}
}
void initializeMisHelper(const Path &path, MisState *states, const Scene* scene,
	size_t nLightPaths, Float misVcWeightFactor, Float misVmWeightFactor, ETransportMode mode){

	for (int i = 0; i < (int)path.vertexCount() - 1; ++i) {
		if (i > 0) states[i] = states[i - 1];
		updateMisHelper(i, path, states[i], scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, mode);
	}
}
void PathSampler::gatherLightPaths(const bool useVC, const bool useVM,
	const float gatherRadius, const int nsample, ImageBlock* lightImage){
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
		m_emitterSubpath.initialize(m_scene, time, EImportance, m_pool);

		/* Perform random walks from the emitter side */
		m_emitterSubpath.randomWalk(m_scene, m_lightPathSampler, m_emitterDepth,
			m_rrDepth, EImportance, m_pool);

		// emitter states
		MisState emitterState, sensorState;
		Spectrum importanceWeight = Spectrum(1.0f);
		for (int s = 1; s < (int)m_emitterSubpath.vertexCount(); ++s) {
			PathVertex
				*vsPred2 = m_emitterSubpath.vertexOrNull(s - 2),
				*vsPred = m_emitterSubpath.vertex(s - 1),
				*vs = m_emitterSubpath.vertex(s);
			PathEdge
				*esPred = m_emitterSubpath.edgeOrNull(s - 2),
				*es = m_emitterSubpath.edgeOrNull(s - 1);

			// Compute the throughput of emitter subpath till this vertex
			importanceWeight *= vsPred->weight[EImportance] * vsPred->rrWeight * es->weight[EImportance];

			// update mis helper and path type
			updateMisHelper(s - 1, m_emitterSubpath, emitterState, m_scene, m_lightPathNum, misVcWeightFactor, misVmWeightFactor, EImportance);

			// store light paths												
			if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon /* don't save backfaced photons */){
				LightVertex lvertex = LightVertex(vs, vsPred, emitterState, importanceWeight);
				m_lightVertices.push_back(lvertex);
				LightVertexExt lvertexExt = LightVertexExt(vs, vsPred, s);
				m_lightVerticesExt.push_back(lvertexExt);
				LightPathNode lnode(vs->getPosition(), m_lightVertices.size() - 1, s);
				m_lightPathTree.push_back(lnode);
			}

			// connect to camera
			if (vs->measure != EDiscrete && lightImage != NULL && useVC){
				Spectrum value = importanceWeight * vs->sampleDirect(m_scene, m_directSampler,
					&vtPred, &vtEdge, &vt, ERadiance);
				if (value.isZero()) continue;

				Point2 samplePos(0.0f);
				vt.getSamplePosition(vs, samplePos);
				bool watchThread = false;
				if (samplePos.x > 50 && samplePos.x < 51 && samplePos.y > 256 && samplePos.y < 257){
					watchThread = true;
				}

				value *= vs->eval(m_scene, vsPred, &vt, EImportance);
				vs->measure = EArea;
				int interactions = m_maxDepth - s;
				if (value.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, &vt, &vtEdge, interactions))
					continue;
				value *= connectionEdge.evalCached(vs, &vt, PathEdge::ETransmittance | PathEdge::ECosineImp);
				Float weight = Path::miWeightVC(watchThread, m_scene, vsPred, vs, &vt, &vtPred,
					s, 1, m_sampleDirect,
					emitterState[EVCM], emitterState[EVC],
					sensorState[EVCM], sensorState[EVC],
					misVmWeightFactor, m_lightPathNum);
				value *= weight;
				if (value.isZero()) continue;

				if (watchThread){
					float fuck = 1.f;
				}

				//value *= invLightPathNum;				
				lightImage->put(samplePos, &value[0]);
			}
		}
		m_lightPathEnds.push_back(m_lightVertices.size());
	}

	// build kdtree
	m_lightPathTree.build(true);

	/* Release any used edges and vertices back to the memory pool */
	m_emitterSubpath.release(m_pool);
}
struct VertexMergingQuery {
	VertexMergingQuery(const Scene *_scene,
		const PathVertex *_vt, const PathVertex *_vtPred, const Vector _wi, const Vector _wiPred,
		const Spectrum &_radianceWeight, const std::vector<LightVertex> &lightVertices,
		int _t, int _maxDepth, Float _misVcWeightFactor, Float _vmNormalization,
		const MisState &_sensorState)
		: scene(_scene), radianceWeight(_radianceWeight), t(_t), maxDepth(_maxDepth),
		vt(_vt), vtPred(_vtPred), wiPred(_wiPred), wi(_wi),
		misVcWeightFactor(_misVcWeightFactor), vmNormalization(_vmNormalization),
		sensorState(_sensorState), m_lightVertices(lightVertices), 
		result(0.f){}

	inline void operator()(const LightPathNode &path) {
		int s = path.data.depth;
		if (maxDepth != -1 && s + t > maxDepth + 2) return;

		size_t vertexIndex = path.data.vertexIndex;
		LightVertex v = m_lightVertices[vertexIndex];
		Vector wo = v.wo;
		MisState emitterState = v.emitterState;
		Spectrum bsdfFactor = vt->eval(scene, wi, wo, ERadiance);
		Spectrum val = v.importanceWeight * radianceWeight * bsdfFactor;

		if (val.isZero()) return;

		Float weightExt = 0.f;
		Float psr2_w = vt->evalPdf(scene, wi, wo, ERadiance, ESolidAngle);
		Float ptr2_w = vt->evalPdf(scene, wo, wi, EImportance, ESolidAngle);
		Float wLight = emitterState[EVCM] * misVcWeightFactor + psr2_w * emitterState[EVM];
		Float wCamera = sensorState[EVCM] * misVcWeightFactor + ptr2_w * sensorState[EVM];
		weightExt = 1.f / (1.f + wLight + wCamera);

		result += val * weightExt * vmNormalization;
	}
	int t, maxDepth;
	Float misVcWeightFactor, vmNormalization;
	MisState sensorState;
	Vector wi, wiPred;	
	const Scene *scene;
	const PathVertex *vt, *vtPred;
	const Spectrum &radianceWeight;
	const std::vector<LightVertex> &m_lightVertices;
	Spectrum result;
};
void PathSampler::sampleCameraPath(const bool useVC, const bool useVM, 
	const float gatherRadius, const Point2i &offset, const size_t cameraPathIndex, SplatList &list) {
	list.clear();

	const Sensor *sensor = m_scene->getSensor();
	size_t nLightPaths = m_lightPathNum;
	const float etaVCM = (M_PI * gatherRadius * gatherRadius) * nLightPaths;
	Float vmNormalization = 1.f / etaVCM;
	Float misVmWeightFactor = useVM ? etaVCM : 0.f;
	Float misVcWeightFactor = useVC ? 1.f / etaVCM : 0.f;	
	Vector2i filmSize = sensor->getFilm()->getSize();

	switch (m_technique) {
	case EBidirectional: {
		/* Uniformly sample a scene time */
		Float time = sensor->getShutterOpen();
		if (sensor->needsTimeSample())
			time = sensor->sampleTime(m_sensorSampler->next1D());

		/* Initialize the path endpoints */
		m_sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);

		if (offset == Point2i(-1))
			m_sensorSubpath.randomWalk(m_scene, m_sensorSampler,
			m_sensorDepth, m_rrDepth, ERadiance, m_pool);
		else
			m_sensorSubpath.randomWalkFromPixel(m_scene, m_sensorSampler,
			m_sensorDepth, offset, m_rrDepth, m_pool);

		/* Compute the combined weights along the two subpaths */
		Spectrum *radianceWeights = (Spectrum *)alloca(m_sensorSubpath.vertexCount()  * sizeof(Spectrum));

		radianceWeights[0] = Spectrum(1.0f);

		for (size_t i = 1; i<m_sensorSubpath.vertexCount(); ++i)
			radianceWeights[i] = radianceWeights[i - 1] *
			m_sensorSubpath.vertex(i - 1)->weight[ERadiance] *
			m_sensorSubpath.vertex(i - 1)->rrWeight *
			m_sensorSubpath.edge(i - 1)->weight[ERadiance];

		bool watchThread = false;
		if (m_sensorSubpath.vertexCount() > 2) {
			Point2 samplePos(0.0f);
			m_sensorSubpath.vertex(1)->getSamplePosition(m_sensorSubpath.vertex(2), samplePos);
			list.append(samplePos, Spectrum(0.0f));

			if (samplePos.x > 50 && samplePos.x < 51 && samplePos.y > 256 && samplePos.y < 257){
				watchThread = true;
			}
		}

		// initialize of MIS helper
		// VCM tech report Eq. 31 - 33		
		MisState *sensorStates = (MisState *)alloca(m_sensorSubpath.vertexCount() * sizeof(MisState));
		initializeMisHelper(m_sensorSubpath, sensorStates, m_scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, ERadiance);

		PathVertex tempEndpoint, tempSample;
		PathEdge tempEdge, connectionEdge;
		Point2 samplePos(0.0f);

		// vertex connection
		if (useVC){
			PathVertex *vs0 = m_pool.allocVertex();
			PathVertex *vsPred0 = m_pool.allocVertex();
			PathEdge *vsEdge0 = m_pool.allocEdge();
			size_t lightPathBegin = (cameraPathIndex == 0) ? 0 : m_lightPathEnds[cameraPathIndex - 1];
			size_t lightPathEnd = m_lightPathEnds[cameraPathIndex];
			for (size_t i = lightPathBegin; i < lightPathEnd + 1; i++){
				int s = 1;
				PathVertex* vsPred = NULL, *vs = NULL;
				PathEdge *vsEdge = NULL;
				Spectrum importanceWeight = Spectrum(1.0f);
				MisState emitterState;

				if (i < lightPathEnd){
					LightVertex lvertex = m_lightVertices[i];
					importanceWeight = lvertex.importanceWeight;
					emitterState = lvertex.emitterState;
					LightVertexExt lvertexExt = m_lightVerticesExt[i];
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
// 				else{
// 					vsEdge = vsEdge0;
// 					vsPred = vsPred0;
// 					vs = vs0;
// 					Spectrum throughput = Spectrum(1.0f);
// 					vsPred->makeEndpoint(m_scene, time, EImportance);
// 					vsPred->sampleNext(m_scene, m_emitterSampler, NULL, NULL, vsEdge,
// 						vs, EImportance, false, &throughput);
// 				}

				int minT = 2, maxT = (int)m_sensorSubpath.vertexCount() - 1;
				if (m_maxDepth != -1)
					maxT = std::min(maxT, m_maxDepth + 1 - s);
				for (int t = minT; t <= maxT; ++t) {
					PathVertex
						*vtPred = m_sensorSubpath.vertexOrNull(t - 1),
						*vt = m_sensorSubpath.vertex(t);
					PathEdge
						*vtEdge = m_sensorSubpath.edgeOrNull(t - 1);

					/* Will receive the path weight of the (s, t)-connection */
					Spectrum value = Spectrum(0.0f);

					if (s == 1){
						value = radianceWeights[t] * vt->sampleDirect(m_scene, m_directSampler,
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
					int remaining = m_maxDepth - depth;					

					// backup original path vertex measure
					uint8_t vsMeasure0 = vs->measure, vtMeasure0 = vt->measure;

					/* Account for the terms of the measurement contribution
					function that are coupled to the connection endpoints */
					if (vt->isSensorSupernode()) {
						/* If possible, convert 'vs' into an sensor sample */
						if (!vs->cast(m_scene, PathVertex::ESensorSample) || vs->isDegenerate())
							continue;

						/* Make note of the changed pixel sample position */
						if (!vs->getSamplePosition(vsPred, samplePos))
							continue;

						value = importanceWeight *
							vs->eval(m_scene, vsPred, vt, EImportance) *
							vt->eval(m_scene, vtPred, vs, ERadiance);
					}
					else if (m_sampleDirect && ((t == 1 && s > 1) || (s == 1 && t > 1))) {
						if (s == 1) {
							if (vt->isDegenerate())
								continue;
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
							value = importanceWeight * vs->sampleDirect(m_scene, m_directSampler,
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

					if (m_excludeDirectIllum && depth <= 2)
						continue;

					/* Account for the terms of the measurement contribution
					function that are coupled to the connection edge */
					if (!sampleDirect)
						value *= connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
					else
						value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
						(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

					MisState sensorState = sensorStates[t - 1];
					Float weight = Path::miWeightVC(watchThread, m_scene, vsPred, vs, vt, vtPred,
						s, t, m_sampleDirect,
						emitterState[EVCM], emitterState[EVC],
						sensorState[EVCM], sensorState[EVC],
						misVmWeightFactor, nLightPaths);
					value *= weight;

					if (weight == 0.0f) continue;

					if (watchThread){
						float fuck = 1.f;
					}

					if (t < 2) {
						list.append(samplePos, value);
					}
					else {
						BDAssert(m_sensorSubpath.vertexCount() > 2);
						list.accum(0, value);
					}
				}
			}
			m_pool.release(vs0);
			m_pool.release(vsPred0);
			m_pool.release(vsEdge0);
		}

		// massive vertex merging
		if (useVM){
			int minT = 2, maxT = (int)m_sensorSubpath.vertexCount() - 1;
			if (m_maxDepth != -1)
				maxT = std::min(maxT, m_maxDepth - 1);

			for (int t = minT; t <= maxT; ++t) {
				PathVertex
					*vtPred2 = m_sensorSubpath.vertexOrNull(t - 2),
					*vtPred = m_sensorSubpath.vertexOrNull(t - 1),
					*vt = m_sensorSubpath.vertex(t);
				PathEdge
					*vtEdge = m_sensorSubpath.edgeOrNull(t - 1);

				if (!vt->isDegenerate()){
					BDAssert(vt->type == PathVertex::ESurfaceInteraction);

					Vector wi = normalize(vtPred->getPosition() - vt->getPosition());
					Vector wiPred = (t == 2) ? Vector(-1.f, -1.f, -1.f) : normalize(vtPred2->getPosition() - vtPred->getPosition());
					VertexMergingQuery query(m_scene, vt, vtPred, wi, wiPred, radianceWeights[t], m_lightVertices,
						t, m_maxDepth, misVcWeightFactor, vmNormalization, sensorStates[t - 1]);

					m_lightPathTree.executeQuery(vt->getPosition(), gatherRadius, query);

					if (!query.result.isZero())
						list.accum(0, query.result);
				}
			}
		}

		/* Release any used edges and vertices back to the memory pool */
		m_sensorSubpath.release(m_pool);
	}
		break;

	default:
		Log(EError, "PathSampler::sample(): invalid technique!");
	}
}



int PathSampler::getConnectionFlag(bool connectionImportance, bool connectionRadiance, bool connectionVisibility,
	bool connectionMIS, bool connectionBSDFs, bool connectionGeometry, bool connectionFull){
	int ret = 0;
	if (connectionImportance) ret |= EConnectImportance;
	if (connectionRadiance) ret |= EConnectRadiance;
	if (connectionVisibility) ret |= EConnectVisibility;
	if (connectionMIS) ret |= EConnectMis;
	if (connectionBSDFs) ret |= EConnectBRDF;
	if (connectionGeometry) ret |= EConnectGeometry;
	if (connectionFull) ret |= EConnectAll;
	return ret;
}

std::string SplatList::toString() const {
	std::ostringstream oss;
	oss << "SplatList[" << endl
		<< "  luminance = " << luminance << "," << endl
		<< "  splats = {" << endl;
	for (size_t i=0; i<splats.size(); ++i) {
			oss << "      " << splats[i].first.toString()
				<< " => " << splats[i].second.toString();
		if (i+1 < splats.size())
			oss << ",";
		oss << endl;
	}
	oss << "  }" << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(PathSampler, false, Object)
MTS_IMPLEMENT_CLASS(SeedWorkUnit, false, WorkUnit)
MTS_NAMESPACE_END
