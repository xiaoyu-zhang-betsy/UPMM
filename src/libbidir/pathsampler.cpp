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

#include <mitsuba/core/sfcurve.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

StatsCounter avgGatherPoints("Unbiased photon mapping", "Avg. number of gather points", EAverage);
StatsCounter maxGatherPoints("Unbiased photon mapping", "Max. number of gather points", EMaximumValue);
StatsCounter numZeroGatherPoints("Unbiased photon mapping", "Percentage of zero gather points", EPercentage);

StatsCounter maxInvpShoots("Unbiased photon mapping", "Max. number of 1/p shoots", EMaximumValue);
StatsCounter numInvpShoots("Unbiased photon mapping", "Total number of 1/p shoots");
StatsCounter avgInvpShoots("Unbiased photon mapping", "Avg. number of 1/p shoots", EAverage);

StatsCounter numClampShoots("Unbiased photon mapping", "Percentage of clamped invp evaluations(bias)", EPercentage);

PathSampler::PathSampler(ETechnique technique, const Scene *scene, Sampler *sensorSampler,
		Sampler *emitterSampler, Sampler *directSampler, int maxDepth, int rrDepth,
		bool excludeDirectIllum,  bool sampleDirect, bool lightImage,
		Sampler *lightPathSampler)
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
						value *= Path::miWeight(m_scene, m_emitterSubpath, &connectionEdge,
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
			value *= Path::miWeight(m_scene, m_emitterSubpath, &connectionEdge,
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
				Float misWeight = Path::miWeight(m_scene, m_emitterSubpath, &connectionEdge,
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
inline Float MisHeuristic(Float pdf) {
	return pdf * pdf;
}
Float miWeightVC(const Scene *scene,
	const PathVertex *vsPred, const PathVertex *vs,
	const PathVertex *vt, const PathVertex *vtPred,
	int s, int t,
	float emitterdVCM, float emitterdVC,
	float sensordVCM, float sensordVC,
	Float misVmWeightFactor, bool isUPM = false){

	Float weight = 0.0f;

	if (s > 1 && t > 1){
		Float psr2_w = vs->evalPdf(scene, vt, vsPred, ERadiance, ESolidAngle);
		Float ptr2_w = vt->evalPdf(scene, vs, vtPred, EImportance, ESolidAngle);
		Float psr1 = vt->evalPdf(scene, vtPred, vs, ERadiance, EArea);
		Float ptr1 = vs->evalPdf(scene, vsPred, vt, EImportance, EArea);

		//Float wCamera = MisHeuristic(ptr1) * ((t == 2 ? 0.f : misVmWeightFactor) + sensordVCM + MisHeuristic(ptr2_w) * sensordVC);
		Float wCamera = MisHeuristic(ptr1) * (misVmWeightFactor + sensordVCM + MisHeuristic(ptr2_w) * sensordVC);
		Float wLight = MisHeuristic(psr1) * ((s == 2 ? 0.f: misVmWeightFactor) + emitterdVCM + MisHeuristic(psr2_w) * emitterdVC);
		//Float wLight = MisHeuristic(psr1) * (misVmWeightFactor + emitterdVCM + MisHeuristic(psr2_w) * emitterdVC);
		weight = 1.f / (1.f + wLight + wCamera);
	}
	else if (t == 1 && s > 1){		
		Float psr2_w = vs->evalPdf(scene, vt, vsPred, ERadiance, ESolidAngle);
		Float psr1 = vt->evalPdf(scene, vtPred, vs, ERadiance, EArea);

		Float wLight;
		if (isUPM){
			wLight = MisHeuristic(psr1) * (((s == 2) ? 0.f : misVmWeightFactor) + emitterdVCM + MisHeuristic(psr2_w) * emitterdVC); // exclude (2,1) path in upm
			//wLight = MisHeuristic(psr1) * (emitterdVCM + MisHeuristic(psr2_w) * emitterdVC); // exclude (2,1) path in upm
		}
		else
			wLight = MisHeuristic(psr1) * (misVmWeightFactor + emitterdVCM + MisHeuristic(psr2_w) * emitterdVC);
		weight = 1.f / (1.f + wLight);
	}
	else if (s == 1 && t > 1){		
		const PositionSamplingRecord &pRec = vs->getPositionSamplingRecord();
		const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
		EMeasure measure = emitter->getDirectMeasure();

		Point p = pRec.p;
		Point p2 = vs->getPosition();

// 		Vector d = normalize(vt->getPosition() - vs->getPosition());
// 		Float pconnect2 = vt->evalPdfDirect(scene, vs, EImportance, measure);		
// 		Float ptrace2 = vsPred->pdf[EImportance];
// 		if (vs->getAbstractEmitter()->needsPositionSample() && vs->getAbstractEmitter()->needsDirectionSample()){
// 			pconnect2 *= absDot(d, vs->getGeometricNormal());
// 		}
// 
// 		if (vs->getAbstractEmitter()->needsDirectionSample()){
// 			DirectionSamplingRecord dRec(d);
// 			ptrace2 *= emitter->pdfDirection(dRec, pRec) * absDot(d, vt->getGeometricNormal());
// 			if (!vs->getAbstractEmitter()->needsPositionSample()){
// 				Float dis = (vt->getPosition() - vs->getPosition()).length();
// 				ptrace2 /= dis * dis;
// 			}
// 		}
// 		Float porig = ptrace2 / pconnect2;

		EMeasure measureTrace = EArea;
		if (emitter->needsPositionSample() && !emitter->needsDirectionSample()){ // directional light
			measureTrace = EDiscrete;
		}		

		Float ptr2_w = vt->evalPdf(scene, vs, vtPred, EImportance, ESolidAngle);
		Float psr1 = vt->evalPdf(scene, vtPred, vs, ERadiance, EArea);

		Float pconnect = vt->evalPdfDirect(scene, vs, EImportance, measure == ESolidAngle ? EArea : measure);
		Float ptrace = vsPred->pdf[EImportance];
		Float ptr1 = vs->evalPdf(scene, vsPred, vt, EImportance, measureTrace); // measure == ESolidAngle ? EArea : measure);
		Float pnew = ptrace / pconnect * ptr1;

		Float wLight = (measure == EDiscrete) ? 0.f : MisHeuristic(psr1 / pconnect);
		Float wCamera;
		if (isUPM){
			//wCamera = MisHeuristic(ptrace / pconnect * ptr1) * ((t == 2 ? 0.f : misVmWeightFactor) + sensordVCM + MisHeuristic(ptr2_w) * sensordVC);
			wCamera = MisHeuristic(ptrace / pconnect * ptr1) * (sensordVCM + MisHeuristic(ptr2_w) * sensordVC);
		}
		else
			wCamera = MisHeuristic(ptrace / pconnect * ptr1) * (misVmWeightFactor + sensordVCM + MisHeuristic(ptr2_w) * sensordVC); // exclude (1,t) path in upm
		weight = 1.f / (1.f + wLight + wCamera);
	}
	else if (s == 1 && t == 1){
		const PositionSamplingRecord &pRec = vs->getPositionSamplingRecord();
		const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
		EMeasure measure = emitter->getDirectMeasure();
		Float pconnect = vt->evalPdfDirect(scene, vs, EImportance, measure == ESolidAngle ? EArea : measure);
		Float psr1 = vt->evalPdf(scene, vtPred, vs, ERadiance, EArea);
		Float wLight = (measure == EDiscrete) ? 0.f : MisHeuristic(psr1 / pconnect);
		weight = 1.f / (1.f + wLight);
	}
	else if (s == 0){
		const PositionSamplingRecord &pRec = vt->getPositionSamplingRecord();
		const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
// 		Float pt_1 = emitter->pdfPosition(pRec);
// 		Vector d = normalize(vtPred->getPosition() - vt->getPosition());
// 		DirectionSamplingRecord dRec(d);
// 		Float pt_2_sigma = emitter->pdfDirection(dRec, pRec);
// 		Float wCamera = MisHeuristic(pt_1) * sensordVCM + MisHeuristic(pt_1 * pt_2_sigma) * sensordVC;

		PathVertex vsTemp;
		vsTemp.makeEndpoint(scene, 0, EImportance);
		EMeasure measure = vt->getAbstractEmitter()->getDirectMeasure();
		Float pconnect = vtPred->evalPdfDirect(scene, vt, EImportance, measure == ESolidAngle ? EArea : measure);
		Float ptrace = vsTemp.evalPdf(scene, NULL, vt, EImportance, measure == ESolidAngle ? EArea : measure);
		Float ptr1_w = vt->evalPdf(scene, &vsTemp, vtPred, EImportance, ESolidAngle);

		//if (t == 3 && misVmWeightFactor > 0.f){
		if (t >= 3 && misVmWeightFactor > 0.f){
			Vector edir = normalize(vt->getPosition() - vtPred->getPosition());
			Float elen = (vt->getPosition() - vtPred->getPosition()).length();
			Float giIn = std::abs(vtPred->isOnSurface() ? dot(edir, vtPred->getGeometricNormal()) : 1) / (elen * elen);
			Float invpi = 1.f / vtPred->pdf[ERadiance];
			sensordVC -= MisHeuristic(giIn * invpi) * misVmWeightFactor;
		}

		Float wCamera = MisHeuristic(pconnect) * sensordVCM + MisHeuristic(ptrace * ptr1_w) * sensordVC;
		
		weight = 1.f / (1.f + wCamera);
	}
	else{
		SLog(EError, "Not implemented miweightVC for s = %d, t = %d", s, t);
	}
	return weight;
}
Float miWeightVM(const Scene *scene, int s, int t,
	const PathVertex *vsPred2, const PathVertex *vsPred, const PathVertex *vs,
	const PathVertex *vt, const PathVertex *vtPred, const PathVertex *vtPred2,
	MisState emitterState, MisState sensorState,	
	Float misVcWeightFactor, bool cameraDirConnection){	

	Float miWeight;
	if (cameraDirConnection){
		Float pt = vtPred->evalPdf(scene, vtPred2, vs, ERadiance, EArea);				
		//Float ps = vsPred->evalPdf(scene, vsPred2, vs, EImportance, EArea);
		Float ps = vsPred->pdf[EImportance];
		Float invpt = (pt > 0.f) ? 1.f / pt : 0.f;
		Float invps = (ps > 0.f) ? 1.f / ps : 0.f;
		if (!_finite(invpt) || invpt == 0.f || !_finite(invps) || invps == 0.f)
			return 0.f;

		Float dVCM = 1.f;
		if (t == 2){
			EMeasure measure = vtPred->getAbstractEmitter()->getDirectMeasure();
			Float pconnect = vs->evalPdfDirect(scene, vtPred, ERadiance, measure == ESolidAngle ? EArea : measure);
			Float ptrace = vtPred2->evalPdf(scene, NULL, vtPred, ERadiance, measure == ESolidAngle ? EArea : measure);
			dVCM = pconnect / ptrace;
		}

		Vector wi = normalize(vtPred->getPosition() - vs->getPosition());
		Vector wo = normalize(vsPred->getPosition() - vs->getPosition());		
		Float psr1_w = vs->evalPdf(scene, wi, wo, ERadiance, ESolidAngle);
		Float ptr1_w = vs->evalPdf(scene, wo, wi, EImportance, ESolidAngle);		
		Float psr2_w = vsPred->evalPdf(scene, vs, vsPred2, ERadiance, ESolidAngle);
		Float ptr2_w = vtPred->evalPdf(scene, vs, vtPred2, EImportance, ESolidAngle);

		Float wLight = emitterState[EVCM] * misVcWeightFactor + MisHeuristic(psr1_w * invps) * (emitterState[EVM] + MisHeuristic(psr2_w) * emitterState[EVMB]);
		Float wCamera = MisHeuristic(dVCM * invpt) * misVcWeightFactor + MisHeuristic(ptr1_w * invpt) * (sensorState[EVM] + MisHeuristic(ptr2_w) * sensorState[EVMB]);
		miWeight = 1.f / (1.f + wLight + wCamera);

		if (!_finite(miWeight) || _isnan(miWeight)){
			SLog(EWarn, "Invalid MIS weight %f, wLight = %f, wCamera = %f, emitterState=[%f, %f, %f, %f], sensorState=[%f, %f, %f, %f], \
					invps = %f, invpt = %f, dVCM = %f, psr1_w = %f, psr2_w = %f, ptr1_w = %f, ptr2_w = %f",
					miWeight, wLight, wCamera, emitterState[EVCM], emitterState[EVC], emitterState[EVM], emitterState[EVMB],
					sensorState[EVCM], sensorState[EVC], sensorState[EVM], sensorState[EVMB],
					invps, invpt, dVCM, psr1_w, psr2_w, ptr1_w, ptr2_w);
		}
	}
	else{
		Float ps = vsPred->evalPdf(scene, vsPred2, vt, EImportance, EArea);
		Float pt = vtPred->pdf[ERadiance]; // evalPdf(scene, vtPred2, vt, ERadiance, EArea);
		Float invps = (ps > 0.f) ? 1.f / ps : 0.f;
		Float invpt = (pt > 0.f) ? 1.f / pt : 0.f;
		if (!_finite(invps) || invps == 0.f || !_finite(invpt) || invpt == 0.f)
			return 0.f;

		Float dVCM = 1.f;
		if (s == 2){
			EMeasure measure = vsPred->getAbstractEmitter()->getDirectMeasure();
			Float pconnect = vt->evalPdfDirect(scene, vsPred, EImportance, measure == ESolidAngle ? EArea : measure);
			Float ptrace = vsPred2->evalPdf(scene, NULL, vsPred, EImportance, measure == ESolidAngle ? EArea : measure);
			dVCM = pconnect / ptrace;
		}

		Vector wi = normalize(vtPred->getPosition() - vt->getPosition());
		Vector wo = normalize(vsPred->getPosition() - vt->getPosition());		
		Float psr1_w = vt->evalPdf(scene, wi, wo, ERadiance, ESolidAngle);
		Float ptr1_w = vt->evalPdf(scene, wo, wi, EImportance, ESolidAngle);
		Float psr2_w = vsPred->evalPdf(scene, vt, vsPred2, ERadiance, ESolidAngle);
		Float ptr2_w = vtPred->evalPdf(scene, vt, vtPred2, EImportance, ESolidAngle);		

		Float wLight = MisHeuristic(dVCM * invps) * misVcWeightFactor + MisHeuristic(psr1_w * invps) * (emitterState[EVM] + MisHeuristic(psr2_w) * emitterState[EVMB]);
		Float wCamera = sensorState[EVCM] * misVcWeightFactor + MisHeuristic(ptr1_w * invpt) * (sensorState[EVM] + MisHeuristic(ptr2_w) * sensorState[EVMB]);
		miWeight = 1.f / (1.f + wLight + wCamera);

		if (!_finite(miWeight) || _isnan(miWeight)){
			SLog(EWarn, "Invalid MIS weight %f, wLight = %f, wCamera = %f, emitterState=[%f, %f, %f, %f], sensorState=[%f, %f, %f, %f], \
					invps = %f, invpt = %f, dVCM = %f, psr1_w = %f, psr2_w = %f, ptr1_w = %f, ptr2_w = %f",
					miWeight, wLight, wCamera, emitterState[EVCM], emitterState[EVC], emitterState[EVM], emitterState[EVMB],
					sensorState[EVCM], sensorState[EVC], sensorState[EVM], sensorState[EVMB],
					invps, invpt, dVCM, psr1_w, psr2_w, ptr1_w, ptr2_w);
		}
	}	

	return miWeight;
}
void updateMisHelper(int i, const Path &path, MisState &state, const Scene* scene,
	size_t nLightPaths, Float misVcWeightFactor, Float misVmWeightFactor, ETransportMode mode, bool isUPM = false){
	if (i == 0){
		state[EVCM] = 0.f;
		state[EVC] = 0.f;
		state[EVM] = 0.f;
		state[EVMB] = 0.f;
	}
	else if (i == 1){
		state[EVCM] = 0.f;
		state[EVC] = 0.f;
		state[EVM] = 0.f;
		state[EVMB] = 0.f;
		PathVertex *v = path.vertex(1);
		PathEdge *e = path.edge(1);
		PathVertex *vNext = path.vertex(2);
		EMeasure measure = v->getAbstractEmitter()->getDirectMeasure();
		Float pconnect = vNext->evalPdfDirect(scene, v, mode, measure == ESolidAngle ? EArea : measure);
		Float ptrace = path.vertex(0)->pdf[mode] * path.edge(0)->pdf[mode];		
		Float p1 = v->pdf[mode];

		if (mode == ERadiance){
			state[EVCM] = MisHeuristic(pconnect / (ptrace * p1));
			state[EVC] = 0.f;
			state[EVM] = 0.f;
		}
		else{
			state[EVCM] = MisHeuristic(pconnect / (ptrace * p1));
			if (measure != EDiscrete){
				Float geoTerm0 = v->isOnSurface() ? std::abs(dot(e->d, v->getGeometricNormal()) / (e->length * e->length)) : 1;
				state[EVC] = MisHeuristic(geoTerm0 / (ptrace * p1));
			}
			else{
				state[EVC] = 0.f;
			}
			if (!isUPM){
				state[EVM] = state[EVC] * misVcWeightFactor;
			}
			else{
				state[EVM] = state[EVC] * misVcWeightFactor * MisHeuristic(p1);
			}
		}
	}
	else{
		PathVertex *v = path.vertex(i);
		PathEdge   *e = path.edge(i);
		if (v->measure == EDiscrete){
			PathVertex *vNext = path.vertex(i + 1);
			state[EVCM] = 0.f;
			state[EVC] *= MisHeuristic(std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / std::abs(vNext->isOnSurface() ? dot(e->d, vNext->getGeometricNormal()) : 1));
			if (!isUPM)
				state[EVM] *= MisHeuristic(std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / std::abs(vNext->isOnSurface() ? dot(e->d, vNext->getGeometricNormal()) : 1));
			else{
				state[EVM] *= MisHeuristic(std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1));
				state[EVMB] *= MisHeuristic(std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1));
			}
		}
		else{
			BDAssert(v->measure == EArea);
			PathVertex *vPred = path.vertex(i - 1);
			PathEdge * ePred = path.edge(i - 1);
			Float giIn = std::abs(v->isOnSurface() ? dot(e->d, v->getGeometricNormal()) : 1) / (e->length * e->length);
			Float giInPred = std::abs(vPred->isOnSurface() ? dot(ePred->d, vPred->getGeometricNormal()) : 1) / (ePred->length * ePred->length);
			Float pi = v->pdf[mode];
			Float pir2_w = (giInPred == 0.f) ? 0.f : (v->pdf[1 - mode] / giInPred);
			Float invpi = 1.f / pi;
			if (isUPM){
				state[EVC] = MisHeuristic(giIn * invpi) * (state[EVCM] + MisHeuristic(pir2_w) * state[EVC]);
				Float piPred = vPred->pdf[mode];
				Float pir2Pred_w = 0.f;
				if (i > 2){
					PathVertex *vPred2 = path.vertex(i - 2);
					PathEdge * ePred2 = path.edge(i - 2);
					Float giInPred2 = std::abs(vPred2->isOnSurface() ? dot(ePred2->d, vPred2->getGeometricNormal()) : 1) / (ePred2->length * ePred2->length);
					pir2Pred_w = (giInPred2 == 0.f) ? 0.f : (vPred->pdf[1 - mode] * ePred2->pdf[1 - mode] / giInPred2);
				}
				state[EVMB] = MisHeuristic(giIn / piPred) * (state[EVM] + MisHeuristic(pir2Pred_w) * state[EVMB]);//state[EVM]
				state[EVM] = MisHeuristic(giIn) * (state[EVCM] * misVcWeightFactor);		

				if (i > 2 || mode == ERadiance){
					state[EVC] += MisHeuristic(giIn * invpi) * misVmWeightFactor;
					state[EVM] += MisHeuristic(giIn);
				}				
				state[EVCM] = MisHeuristic(invpi);
			}
			else{				
				state[EVC] = MisHeuristic(giIn * invpi) * (state[EVCM] + misVmWeightFactor + MisHeuristic(pir2_w) * state[EVC]);
				state[EVM] = MisHeuristic(giIn * invpi) * (1.f + state[EVCM] * misVcWeightFactor + MisHeuristic(pir2_w) * state[EVM]);
				state[EVCM] = MisHeuristic(invpi);
			}			

// 			if (i == 2 && isUPM){
// 				state[EVC] = MisHeuristic(giIn) * (state[EVCM] + MisHeuristic(pir2_w) * state[EVC]);
// 				state[EVM] = MisHeuristic(giIn) * (state[EVCM] * misVcWeightFactor + MisHeuristic(pir2_w) * state[EVM]);
// 			}
// 			else{
// 				state[EVC] = MisHeuristic(giIn) * (state[EVCM] + misVmWeightFactor + MisHeuristic(pir2_w) * state[EVC]);
// 				state[EVM] = MisHeuristic(giIn) * (1.f + state[EVCM] * misVcWeightFactor + MisHeuristic(pir2_w) * state[EVM]);
// 			}
// 			state[EVCM] = MisHeuristic(1 / pi);
// 			state[EVC] *= state[EVCM];
// 			if (mode == EImportance || !isUPM){
// 				state[EVM] *= state[EVCM];
// 			}
		}
	}
}
void initializeMisHelper(const Path &path, MisState *states, const Scene* scene,
	size_t nLightPaths, Float misVcWeightFactor, Float misVmWeightFactor, ETransportMode mode, bool isUPM = false){

	for (int i = 0; i < (int)path.vertexCount() - 1; ++i) {
		if (i > 0) states[i] = states[i - 1];
		updateMisHelper(i, path, states[i], scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, mode, isUPM);
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
	m_lightPathNum = nsample;
	Float etaVCM = (M_PI * gatherRadius * gatherRadius) * m_lightPathNum;
	Float invLightPathNum = 1.f / m_lightPathNum;
	Float misVmWeightFactor = useVM ? MisHeuristic(etaVCM) : 0.f;
	Float misVcWeightFactor = useVC ? MisHeuristic(1.f / etaVCM) : 0.f;
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
			LightVertex lvertex = LightVertex(vs, vsPred, emitterState, importanceWeight);
			m_lightVertices.push_back(lvertex);
			LightVertexExt lvertexExt = LightVertexExt(vs, vsPred, s);
			m_lightVerticesExt.push_back(lvertexExt);
			if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon){
				LightPathNode lnode(vs->getPosition(), m_lightVertices.size() - 1, s);
				m_lightPathTree.push_back(lnode);
			}

			// connect to camera
			if (vs->measure != EDiscrete && lightImage != NULL && useVC){
				Spectrum value = importanceWeight * vs->sampleDirect(m_scene, m_directSampler,
					&vtPred, &vtEdge, &vt, ERadiance);
				if (value.isZero()) continue;								

				value *= vs->eval(m_scene, vsPred, &vt, EImportance);
				vs->measure = EArea;
				int interactions = m_maxDepth - s;
				if (value.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, &vt, &vtEdge, interactions))
					continue;
				value *= connectionEdge.evalCached(vs, &vt, PathEdge::ETransmittance | PathEdge::ECosineImp);
				Float weight = miWeightVC(m_scene, vsPred, vs, &vt, &vtPred,
					s, 1,
					emitterState[EVCM], emitterState[EVC],
					sensorState[EVCM], sensorState[EVC],
					misVmWeightFactor);
				value *= weight;
				if (value.isZero()) continue;
				
				Point2 samplePos(0.0f);
				vt.getSamplePosition(vs, samplePos);
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
		Float wLight = emitterState[EVCM] * misVcWeightFactor + MisHeuristic(psr2_w) * emitterState[EVM];
		Float wCamera = sensorState[EVCM] * misVcWeightFactor + MisHeuristic(ptr2_w) * sensorState[EVM];
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
void PathSampler::sampleSplatsVCM(const bool useVC, const bool useVM, 
	const float gatherRadius, const Point2i &offset, const size_t cameraPathIndex, SplatList &list) {
	list.clear();

	const Sensor *sensor = m_scene->getSensor();
	size_t nLightPaths = m_lightPathNum;
	const float etaVCM = (M_PI * gatherRadius * gatherRadius) * nLightPaths;
	Float vmNormalization = 1.f / etaVCM;
	Float misVmWeightFactor = useVM ? MisHeuristic(etaVCM) : 0.f;
	Float misVcWeightFactor = useVC ? MisHeuristic(1.f / etaVCM) : 0.f;

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

		Point2 initialSamplePos(0.0f);
		if (m_sensorSubpath.vertexCount() > 2) {
			Point2 samplePos(0.0f);
			m_sensorSubpath.vertex(1)->getSamplePosition(m_sensorSubpath.vertex(2), samplePos);
			list.append(samplePos, Spectrum(0.0f));
			initialSamplePos = samplePos;
		}

		// initialize of MIS helper
		// VCM tech report Eq. 31 - 33		
		MisState *sensorStates = (MisState *)alloca(m_sensorSubpath.vertexCount() * sizeof(MisState));
		initializeMisHelper(m_sensorSubpath, sensorStates, m_scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, ERadiance);

		// massive vertex merging
		if (useVM){
			int minT = 2, maxT = (int)m_sensorSubpath.vertexCount() - 1;
			if (m_maxDepth != -1)
				maxT = std::min(maxT, m_maxDepth);

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

		PathVertex tempEndpoint, tempSample;
		PathEdge tempEdge, connectionEdge;		
		// vertex connection
		if (useVC){
			PathVertex *vs0 = m_pool.allocVertex();
			PathVertex *vsPred0 = m_pool.allocVertex();
			size_t lightPathBegin = (cameraPathIndex == 0) ? 0 : m_lightPathEnds[cameraPathIndex - 1];
			size_t lightPathEnd = m_lightPathEnds[cameraPathIndex];
			Point2 samplePos(0.0f);
			PathEdge connectionEdge;
			for (size_t i = lightPathBegin; i < lightPathEnd + 2; i++){
				int s = lightPathEnd + 1 - i;
				PathVertex* vsPred = vsPred0, *vs = vs0;
				Spectrum importanceWeight = Spectrum(1.0f);
				MisState emitterState;

				memset(vsPred, 0, sizeof(PathVertex));
				memset(vs, 0, sizeof(PathVertex));
				if (i < lightPathEnd){
					LightVertex lvertex = m_lightVertices[i];
					importanceWeight = lvertex.importanceWeight;
					emitterState = lvertex.emitterState;
					s = m_lightVerticesExt[i].depth;
					if (s == 1) continue;

					m_lightVerticesExt[i].expand(vs);
					m_lightVerticesExt[i - 1].expand(vsPred);
				}

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

					RestoreMeasureHelper rmh1(vt);

					/* Will be set to true if direct sampling was used */
					bool sampleDirect = false;

					/* Number of edges of the combined subpaths */
					int depth = s + t - 1;

					/* Allowed remaining number of ENull vertices that can
					be bridged via pathConnect (negative=arbitrarily many) */
					int remaining = m_maxDepth - depth;

					/* Account for the terms of the measurement contribution
					function that are coupled to the connection endpoints */
					if (s == 0){
						/* If possible, convert 'vt' into an emitter sample */
						if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
							continue;
						vs = vs0;
						vs->makeEndpoint(m_scene, time, EImportance);
						value = radianceWeights[t] *
							vs->eval(m_scene, NULL, vt, EImportance) *
							vt->eval(m_scene, vtPred, vs, ERadiance);
					}
					else if (s == 1) {
						if (vt->isDegenerate())
							continue;
						/* Generate a position on an emitter using direct sampling */
						value = radianceWeights[t] * vt->sampleDirect(m_scene, m_directSampler,
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
					if (vt->isSensorSample())
						if (!vt->getSamplePosition(vs, samplePos))
							continue;

					/* Account for the terms of the measurement contribution
					function that are coupled to the connection edge */
					if (!sampleDirect)
						value *= connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
					else
						value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
						(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

					MisState sensorState = sensorStates[t - 1];
					Float miWeight = miWeightVC(m_scene, vsPred, vs, vt, vtPred,
						s, t,
						emitterState[EVCM], emitterState[EVC],
						sensorState[EVCM], sensorState[EVC],
						misVmWeightFactor, false);

					value *= miWeight;

					if (value.isZero()) continue;

					if (t < 2) {
						list.append(samplePos, value);
					}
					else {
						list.accum(0, value);
					}
				}
			}
			m_pool.release(vs0);
			m_pool.release(vsPred0);
		}

		/* Release any used edges and vertices back to the memory pool */
		m_sensorSubpath.release(m_pool);
	}
		break;

	default:
		Log(EError, "PathSampler::sample(): invalid technique!");
	}
}

void PathSampler::gatherLightPathsUPM(const bool useVC, const bool useVM,
	const float gatherRadius, const int nsample, UPMWorkResult *wr, Float rejectionProb){
	const Sensor *sensor = m_scene->getSensor();
	m_lightPathTree.clear();
	m_lightVertices.clear();
	m_lightVerticesExt.clear();
	m_lightPathEnds.clear();

	PathVertex vtPred, vt;
	PathEdge vtEdge, connectionEdge;

	Float time = sensor->getShutterOpen();
	m_lightPathNum = nsample;
	Float etaVCM = (M_PI * gatherRadius * gatherRadius) * m_lightPathNum * (1.f - rejectionProb);
	Float invLightPathNum = 1.f / m_lightPathNum;
	Float misVmWeightFactor = useVM ? MisHeuristic(etaVCM) : 0.f;
	Float misVcWeightFactor = useVC ? MisHeuristic(1.f / etaVCM) : 0.f;
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
			updateMisHelper(s - 1, m_emitterSubpath, emitterState, m_scene, m_lightPathNum, misVcWeightFactor, misVmWeightFactor, EImportance, true);

			// store light paths												
			//if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon /* don't save backfaced photons */){
			{
				LightVertex lvertex = LightVertex(vs, vsPred, emitterState, importanceWeight);
				m_lightVertices.push_back(lvertex);
				LightVertexExt lvertexExt = LightVertexExt(vs, vsPred, s);
				m_lightVerticesExt.push_back(lvertexExt);
				if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon){
					LightPathNode lnode(vs->getPosition(), m_lightVertices.size() - 1, s);
					m_lightPathTree.push_back(lnode);
				}
			}

			// connect to camera
			if (vs->measure != EDiscrete && wr != NULL && useVC){
				Point2 samplePos(0.0f);
				Spectrum value = importanceWeight * vs->sampleDirect(m_scene, m_directSampler,
					&vtPred, &vtEdge, &vt, ERadiance);
				if (value.isZero() || !vt.getSamplePosition(vs, samplePos)) continue;

				value *= vs->eval(m_scene, vsPred, &vt, EImportance);
				vs->measure = EArea;
				int interactions = m_maxDepth - s;
				if (value.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, &vt, &vtEdge, interactions))
					continue;
				value *= connectionEdge.evalCached(vs, &vt, PathEdge::ETransmittance | PathEdge::ECosineImp);
				Float weight = miWeightVC(m_scene, vsPred, vs, &vt, &vtPred,
					s, 1,
					emitterState[EVCM], emitterState[EVC],
					sensorState[EVCM], sensorState[EVC],
					misVmWeightFactor, true);

#if UPM_DEBUG == 1
				wr->putDebugSample(s, 1, samplePos, value * weight);
				//wr->putDebugSampleM(s, 1, samplePos, value);
#endif

				value *= weight;
				if (value.isZero()) continue;
				
				vt.getSamplePosition(vs, samplePos);
				wr->putSample(samplePos, &value[0]);
			}
		}
		m_lightPathEnds.push_back(m_lightVertices.size());
	}

	// build kdtree
	m_lightPathTree.build(true);

	/* Release any used edges and vertices back to the memory pool */
	m_emitterSubpath.release(m_pool);
}

void PathSampler::gatherCameraPathsUPM(const bool useVC, const bool useVM, 
		const float gatherRadius){

	const Sensor *sensor = m_scene->getSensor();
	m_cameraPathTree.clear();
	m_cameraVertices.clear();
	m_cameraVerticesExt.clear();
	
	HilbertCurve2D<int> hilbertCurve;
	TVector2<int> filmSize(sensor->getFilm()->getCropSize());
	hilbertCurve.initialize(filmSize);
	
	Float time = sensor->getShutterOpen();
	Float etaVCM = (M_PI * gatherRadius * gatherRadius) * m_lightPathNum;
	Float invLightPathNum = 1.f / m_lightPathNum;
	Float misVmWeightFactor = useVM ? MisHeuristic(etaVCM) : 0.f;
	Float misVcWeightFactor = useVC ? MisHeuristic(1.f / etaVCM) : 0.f;
	for (size_t k = 0; k < hilbertCurve.getPointCount(); ++k) {
		Point2i offset = Point2i(hilbertCurve[k]);
		// m_sampler->generate(offset); // TODO: what's this??

		/* Initialize the path endpoints */
		m_sensorSubpath.initialize(m_scene, time, ERadiance, m_pool);

		/* Perform random walks from the emitter side */
		m_sensorSubpath.randomWalkFromPixel(m_scene, m_lightPathSampler,
			m_sensorDepth, offset, m_rrDepth, m_pool);

		Point2 initialSamplePos(0.0f);
		if (m_sensorSubpath.vertexCount() > 2) {
			Point2 samplePos(0.0f);
			m_sensorSubpath.vertex(1)->getSamplePosition(m_sensorSubpath.vertex(2), samplePos);
			initialSamplePos = samplePos;
		}

		// emitter states
		MisState sensorState;
		Spectrum radianceWeight = Spectrum(1.0f);
		for (int t = 1; t < (int)m_sensorSubpath.vertexCount(); ++t) {
			PathVertex
				*vtPred2 = m_sensorSubpath.vertexOrNull(t - 2),
				*vtPred = m_sensorSubpath.vertex(t - 1),
				*vt = m_sensorSubpath.vertex(t);
			PathEdge
				*etPred = m_sensorSubpath.edgeOrNull(t - 2),
				*et = m_sensorSubpath.edgeOrNull(t - 1);

			// Compute the throughput of emitter subpath till this vertex
			radianceWeight *= vtPred->weight[ERadiance] * vtPred->rrWeight * et->weight[ERadiance];

			// update mis helper and path type
			updateMisHelper(t - 1, m_sensorSubpath, sensorState, m_scene, m_lightPathNum, misVcWeightFactor, misVmWeightFactor, ERadiance, true);

			// store light paths												
			//if (s > 1 && vs->measure != EDiscrete && dot(es->d, -vs->getGeometricNormal()) > Epsilon /* don't save backfaced photons */){
			{
				LightVertex lvertex = LightVertex(vt, vtPred, sensorState, radianceWeight, ERadiance, initialSamplePos);
				m_cameraVertices.push_back(lvertex);
				LightVertexExt lvertexExt = LightVertexExt(vt, vtPred, t);
				m_cameraVerticesExt.push_back(lvertexExt);
				if (t > 1 && vt->measure != EDiscrete && dot(et->d, -vt->getGeometricNormal()) > Epsilon){
					LightPathNode lnode(vt->getPosition(), m_cameraVertices.size() - 1, t);
					m_cameraPathTree.push_back(lnode);
				}
			}			
		}
	}

	// build kdtree
	m_cameraPathTree.build(true);

	/* Release any used edges and vertices back to the memory pool */
	m_sensorSubpath.release(m_pool);
}

void PathSampler::sampleSplatsUPM(UPMWorkResult *wr,
	const float gatherRadius, const Point2i &offset, 
	const size_t cameraPathIndex, SplatList &list, bool useVC, bool useVM, 
	Float rejectionProb, size_t clampThreshold) {
	list.clear();

	const Sensor *sensor = m_scene->getSensor();
	size_t nLightPaths = m_lightPathNum;
	const float etaVCM = (M_PI * gatherRadius * gatherRadius) * nLightPaths * (1.f - rejectionProb);
	Float vmNormalization = 1.f / etaVCM;
	Float misVmWeightFactor = useVM ? MisHeuristic(etaVCM) : 0.f;
	Float misVcWeightFactor = useVC ? MisHeuristic(1.f / etaVCM) : 0.f;

	Float squareRadiusPi = M_PI * gatherRadius * gatherRadius;
	Float invLightPaths = 1.f / (float)nLightPaths;

	PathVertex tempSample, tempEndpoint;
	PathEdge tempEdge;
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
		Point2 initialSamplePos(0.0f);
		if (m_sensorSubpath.vertexCount() > 2) {
			Point2 samplePos(0.0f);
			m_sensorSubpath.vertex(1)->getSamplePosition(m_sensorSubpath.vertex(2), samplePos);
			list.append(samplePos, Spectrum(0.0f));
			initialSamplePos = samplePos;
		}

		// initialize of MIS helper
		// VCM tech report Eq. 31 - 33		
		MisState *sensorStates = (MisState *)alloca(m_sensorSubpath.vertexCount() * sizeof(MisState));
		initializeMisHelper(m_sensorSubpath, sensorStates, m_scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, ERadiance, true);

		// massive vertex merging
		if (useVM){
			PathVertex *vs0 = m_pool.allocVertex();
			PathVertex *vsPred0 = m_pool.allocVertex();
			PathVertex *vsPred20 = m_pool.allocVertex();
			PathVertex *vs = NULL, *vsPred = NULL, *vsPred2 = NULL;
			PathEdge connectionEdge;
			PathVertex *succVertex = m_pool.allocVertex();
			PathEdge *succEdge = m_pool.allocEdge();
			Point2 samplePos(0.0f);
			std::vector<uint32_t> searchResults;
			//std::vector<Point> searchPos;
			std::vector<uint32_t> acceptCnt;
			std::vector<size_t> shootCnt;

			int minT = 2; int minS = 3;

			int maxT = (int)m_sensorSubpath.vertexCount() - 1;
			if (m_maxDepth != -1)
				maxT = std::min(maxT, m_maxDepth);
			for (int t = minT; t <= maxT; ++t) {
				PathVertex
					*vtPred2 = m_sensorSubpath.vertex(t - 2),
					*vtPred = m_sensorSubpath.vertex(t - 1),
					*vt = m_sensorSubpath.vertex(t);
				PathEdge
					*predEdge = m_sensorSubpath.edge(t - 2);

				if (!vt->isDegenerate()){
					BDAssert(vt->type == PathVertex::ESurfaceInteraction);

					searchResults.clear();
					m_lightPathTree.search(vt->getPosition(), gatherRadius, searchResults);
					avgGatherPoints.incrementBase();
					avgGatherPoints += searchResults.size();
					maxGatherPoints.recordMaximum(searchResults.size());
					numZeroGatherPoints.incrementBase();
					if (searchResults.size() == 0){
						++numZeroGatherPoints;
						continue;
					}

					//searchPos.resize(searchResults.size());
					acceptCnt.resize(searchResults.size());
					shootCnt.resize(searchResults.size());
					for (int i = 0; i < searchResults.size(); i++){
						acceptCnt[i] = 0;
						shootCnt[i] = 0;
// 						LightPathNode node = m_lightPathTree[searchResults[i]];
// 						size_t vertexIndex = node.data.vertexIndex;
// 						LightVertexExt lvertexExt = m_lightVerticesExt[vertexIndex];
// 						searchPos[i] = lvertexExt.position;
					}

					// evaluate sampling domain pdf normalization
					int shareShootThreshold = 32;
					Float invBrdfIntegral = 1.f;
					bool shareShoot = false;
// 					if (searchPos.size() > shareShootThreshold){
// 						Vector4 smplBBox = Vector4(0.f);
// 						Vector4 smplBBoxDiff = Vector4(0.f);
// 						Float brdfIntegral = vtPred->gatherAreaPdf(vt->getPosition(), gatherRadius * 2.f, vtPred2, smplBBox, &smplBBoxDiff);
// 						if (brdfIntegral > 0.f){
// 							shareShoot = true;
// 							invBrdfIntegral = 1.f / brdfIntegral;
// 							size_t totalShoot = 0, targetShoot = 1;
// 							uint32_t finishCnt = 0;
// 							Float distSquared = gatherRadius * gatherRadius;
// 							while (finishCnt < searchResults.size() && totalShoot < clampThreshold){
// 								totalShoot++;
// 
// 								// restricted sampling evaluation shoots
// 								if (!vtPred->sampleShoot(m_scene, m_sensorSampler, vtPred2, predEdge, succEdge, succVertex, ERadiance, vt->getPosition(), gatherRadius * 2.f, smplBBox, smplBBoxDiff))
// 									continue;
// 
// 								Point pshoot = succVertex->getPosition();
// 								for (int i = 0; i < searchPos.size(); i++){
// 									if (shootCnt[i] > 0) continue;
// 									Float pointDistSquared = (succVertex->getPosition() - searchPos[i]).lengthSquared();
// 									if (pointDistSquared < distSquared){
// 										acceptCnt[i]++;
// 										if (acceptCnt[i] == targetShoot){
// 											shootCnt[i] = totalShoot;
// 											finishCnt++;
// 										}
// 									}
// 								}
// 							}
// 							avgInvpShoots.incrementBase();
//							avgInvpShoots += totalShoot;
// 							maxInvpShoots.recordMaximum(totalShoot);
// 							numInvpShoots += totalShoot;
// 							numClampShoots.incrementBase(searchResults.size());
// 							if (finishCnt < searchResults.size()){
// 								numClampShoots += searchResults.size() - finishCnt;
// 							}
// 						}
// 					}					

					MisState sensorState = sensorStates[t - 1];
					for (int k = 0; k < searchResults.size(); k++){
						LightPathNode node = m_lightPathTree[searchResults[k]];
						int s = node.data.depth;
						if (m_maxDepth != -1 && s + t > m_maxDepth + 2 || s < minS) continue;

						if (s == 2 && t == 2) continue;

						if (m_sensorSampler->next1D() < rejectionProb) continue;

						size_t vertexIndex = node.data.vertexIndex;
						LightVertex vi = m_lightVertices[vertexIndex];
						LightVertex viPred = m_lightVertices[vertexIndex - 1];
						MisState emitterState = vi.emitterState;

						vs = vs0; vsPred = vsPred0; vsPred2 = vsPred20;
						m_lightVerticesExt[vertexIndex].expand(vs);						
						m_lightVerticesExt[vertexIndex - 1].expand(vsPred);						
						if (s > 2)
							m_lightVerticesExt[vertexIndex - 2].expand(vsPred2);
						else
							vsPred2->makeEndpoint(m_scene, time, EImportance);

						// decide the direction to do connection
						bool cameraDirConnection = true;						
						Float sBandwidth = 0.f, tBandwidth = 0.f;
						if (vtPred->isSensorSample()){
							tBandwidth = 10000.f;
						}
						else{
							BDAssert(vtPred->isSurfaceInteraction());
							const Intersection &its = vtPred->getIntersection();
							const BSDF *bsdf = its.getBSDF();
							tBandwidth = bsdf->getBandwidth();
						}
						if (vsPred->isEmitterSample()){
							PositionSamplingRecord &pRec = vsPred->getPositionSamplingRecord();
							const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
							if (!emitter->needsDirectionSample())
								sBandwidth = 99999.f;
							else
								sBandwidth = 0.f;
						}
						else{
							BDAssert(vsPred->isSurfaceInteraction());
							const Intersection &its = vsPred->getIntersection();
							const BSDF *bsdf = its.getBSDF();
							sBandwidth = bsdf->getBandwidth();
						}
 						if (sBandwidth < tBandwidth)
 							cameraDirConnection = false;						

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
							int interactions = m_maxDepth - s - t + 1;
							if (contrib.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, vtPred, NULL, interactions))
								continue;
							contrib *= connectionEdge.evalCached(vs, vtPred, PathEdge::EGeneralizedGeometricTerm);							
						}
						else{							
							contrib = radianceWeights[t] * viPred.importanceWeight * invLightPaths;
							contrib *= vt->eval(m_scene, vtPred, vsPred, ERadiance) *	vsPred->eval(m_scene, vsPred2, vt, EImportance);
							int interactions = m_maxDepth - s - t + 1;
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
							Vector4 smplBBox = Vector4(0.f);
							Vector4 smplBBoxDiff = Vector4(0.f);

							Float brdfIntegral;
							if (cameraDirConnection)
								brdfIntegral = vtPred->gatherAreaPdf(vs->getPosition(), gatherRadius, vtPred2, smplBBox, &smplBBoxDiff);
							else
								brdfIntegral = vsPred->gatherAreaPdf(vt->getPosition(), gatherRadius, vsPred2, smplBBox, &smplBBoxDiff);

							if (brdfIntegral == 0.f) continue;
							invBrdfIntegral = 1.f / brdfIntegral;
							size_t totalShoot = 0, acceptedShoot = 0, targetShoot = 1;
							Float distSquared = gatherRadius * gatherRadius;
							while (totalShoot < clampThreshold){
								totalShoot++;

								// restricted sampling evaluation shoots
								Float pointDistSquared;
								if (cameraDirConnection){
									if (!vtPred->sampleShoot(m_scene, m_sensorSampler, vtPred2, predEdge, succEdge, succVertex, ERadiance, vs->getPosition(), gatherRadius, smplBBox, smplBBoxDiff))
										continue;
									pointDistSquared = (succVertex->getPosition() - vs->getPosition()).lengthSquared();
								}
								else{
									if (!vsPred->sampleShoot(m_scene, m_emitterSampler, vsPred2, predEdge, succEdge, succVertex, EImportance, vt->getPosition(), gatherRadius, smplBBox, smplBBoxDiff))
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
							avgInvpShoots.incrementBase();
							avgInvpShoots += totalShoot;
							maxInvpShoots.recordMaximum(totalShoot);
							numInvpShoots += totalShoot;

							numClampShoots.incrementBase();
							if (totalShoot >= clampThreshold)
								++numClampShoots;

							invp = (acceptedShoot > 0) ? (Float)(totalShoot) / (Float)(acceptedShoot)* invBrdfIntegral : 0;
							contrib *= invp;
						}

						// accumulate to image
						if (contrib.isZero()) continue;

						// MIS weighting						
						Float miWeight = miWeightVM(m_scene, s, t, vsPred2, vsPred, vs, vt, vtPred, vtPred2,
							emitterState, sensorState, misVcWeightFactor, cameraDirConnection);

#if UPM_DEBUG == 1
						if (cameraDirConnection){
							wr->putDebugSample(s, t - 1, samplePos, contrib * miWeight);
							//wr->putDebugSampleM(s, t - 1, samplePos, contrib);
						}
						else{
							wr->putDebugSample(s - 1, t, samplePos, contrib * miWeight);
							//wr->putDebugSampleM(s - 1, t, samplePos, contrib);
						}
#endif

						contrib *= miWeight;
						
#ifdef UPM_DEBUG_HARD
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
			m_pool.release(vs0);
			m_pool.release(vsPred0);
			m_pool.release(vsPred20);			
			m_pool.release(succVertex);
			m_pool.release(succEdge);
		}

		// vertex connection	 
		if (useVC){
			PathVertex *vs0 = m_pool.allocVertex();
			PathVertex *vsPred0 = m_pool.allocVertex();
			size_t lightPathBegin = (cameraPathIndex == 0) ? 0 : m_lightPathEnds[cameraPathIndex - 1];
			size_t lightPathEnd = m_lightPathEnds[cameraPathIndex];
			Point2 samplePos(0.0f);			
			PathEdge connectionEdge;
			for (size_t i = lightPathBegin; i < lightPathEnd + 2; i++){
				int s = lightPathEnd + 1 - i;
				PathVertex* vsPred = vsPred0, *vs = vs0;
				Spectrum importanceWeight = Spectrum(1.0f);
				MisState emitterState;

				memset(vsPred, 0, sizeof(PathVertex));
				memset(vs, 0, sizeof(PathVertex));
				if (i < lightPathEnd){
					LightVertex lvertex = m_lightVertices[i];
					importanceWeight = lvertex.importanceWeight;
					emitterState = lvertex.emitterState;
					s = m_lightVerticesExt[i].depth;
					if (s == 1) continue;

					m_lightVerticesExt[i].expand(vs);
					m_lightVerticesExt[i - 1].expand(vsPred);
				}

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

					RestoreMeasureHelper rmh1(vt);

					/* Will be set to true if direct sampling was used */
					bool sampleDirect = false;

					/* Number of edges of the combined subpaths */
					int depth = s + t - 1;

					/* Allowed remaining number of ENull vertices that can
					be bridged via pathConnect (negative=arbitrarily many) */
					int remaining = m_maxDepth - depth;

					/* Account for the terms of the measurement contribution
					function that are coupled to the connection endpoints */
					if (s == 0){
						/* If possible, convert 'vt' into an emitter sample */
						if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
							continue;
						vs = vs0;
						vs->makeEndpoint(m_scene, time, EImportance);
						value = radianceWeights[t] *
							vs->eval(m_scene, NULL, vt, EImportance) *
							vt->eval(m_scene, vtPred, vs, ERadiance);
					}else if (s == 1) {
						if (vt->isDegenerate())
							continue;
						/* Generate a position on an emitter using direct sampling */
						value = radianceWeights[t] * vt->sampleDirect(m_scene, m_directSampler,
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
					if (vt->isSensorSample())
						if (!vt->getSamplePosition(vs, samplePos))
							continue;

					/* Account for the terms of the measurement contribution
					function that are coupled to the connection edge */
					if (!sampleDirect)
						value *= connectionEdge.evalCached(vs, vt, PathEdge::EGeneralizedGeometricTerm);
					else
						value *= connectionEdge.evalCached(vs, vt, PathEdge::ETransmittance |
						(s == 1 ? PathEdge::ECosineRad : PathEdge::ECosineImp));

					MisState sensorState = sensorStates[t - 1];
					Float miWeight = miWeightVC(m_scene, vsPred, vs, vt, vtPred,
						s, t,
						emitterState[EVCM], emitterState[EVC],
						sensorState[EVCM], sensorState[EVC],
						misVmWeightFactor, true);

#if UPM_DEBUG == 1
					wr->putDebugSample(s, t, samplePos, value * miWeight);
					//wr->putDebugSampleM(s, t, samplePos, value);
#endif

					value *= miWeight;

					if (value.isZero()) continue;

#ifdef UPM_DEBUG_HARD
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
			m_pool.release(vs0);
			m_pool.release(vsPred0);
		}		

		/* Release any used edges and vertices back to the memory pool */
		m_sensorSubpath.release(m_pool);		
	}
		break;

	default:
		Log(EError, "PathSampler::sample(): invalid technique!");
	}
}

Float PathSampler::generateSeedsExtend(const bool useVC, const bool useVM, const float gatherRadius,
	size_t sampleCount, size_t seedCount, std::vector<PathSeed> &seeds) {
	Log(EInfo, "Integrating luminance values over the image plane ("
		SIZE_T_FMT " samples)..", sampleCount);

	BDAssert(m_sensorSampler == m_emitterSampler);
	BDAssert(m_sensorSampler->getClass()->derivesFrom(MTS_CLASS(ReplayableSampler)));

	ref<Timer> timer = new Timer();
	std::vector<PathSeed> tempSeeds;
	tempSeeds.reserve(sampleCount);

	SplatList splatList;
	Float luminance;
// 	PathCallback callback = boost::bind(&seedCallback,
// 		boost::ref(tempSeeds), NULL, boost::ref(luminance),
// 		_1, _2, _3, _4);

	double mean = 0.0, variance = 0.0;
	for (size_t i = 0; i<sampleCount; ++i) {
		size_t seedIndex = tempSeeds.size();
		size_t sampleIndex = m_sensorSampler->getSampleIndex();
		luminance = 0.0f;

		/* Run the path sampling strategy */
		sampleSplatsExtend(useVC, useVM, gatherRadius, Point2i(-1), splatList);
		luminance = splatList.luminance;
		splatList.normalize(NULL);

		/* Coarse seed granularity (e.g. for PSSMLT) */
		if (luminance != 0)
			tempSeeds.push_back(PathSeed(sampleIndex, luminance));

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

void PathSampler::sampleSplatsExtend(const bool useVC, const bool useVM, const float gatherRadius,
	const Point2i &offset, SplatList &list) {
	const Sensor *sensor = m_scene->getSensor();

	size_t nLightPaths = m_lightPathNum;
	Float invLightPaths = 1.f / (float)nLightPaths;
	const float etaVCM = (M_PI * gatherRadius * gatherRadius) * nLightPaths;
	Float vmNormalization = 1.f / etaVCM;
	Float misVmWeightFactor = useVM ? MisHeuristic(etaVCM) : 0.f;
	Float misVcWeightFactor = useVC ? MisHeuristic(1.f / etaVCM) : 0.f;

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

		Point2 initialSamplePos(0.0f);
		if (m_sensorSubpath.vertexCount() > 2) {
			Point2 samplePos(0.0f);
			m_sensorSubpath.vertex(1)->getSamplePosition(m_sensorSubpath.vertex(2), samplePos);
			list.append(samplePos, Spectrum(0.0f));
			initialSamplePos = samplePos;
		}

		// initialize of MIS helper
		// VCM tech report Eq. 31 - 33		
		MisState *sensorStates = (MisState *)alloca(m_sensorSubpath.vertexCount() * sizeof(MisState));
		initializeMisHelper(m_sensorSubpath, sensorStates, m_scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, ERadiance, true);
		MisState *emitterStates = (MisState *)alloca(m_emitterSubpath.vertexCount() * sizeof(MisState));
		initializeMisHelper(m_emitterSubpath, emitterStates, m_scene, nLightPaths, misVcWeightFactor, misVmWeightFactor, EImportance, true);

		PathVertex tempEndpoint, tempSample;
		PathEdge tempEdge, connectionEdge;
		Point2 samplePos(0.0f);

		if (useVC){
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

					/* Account for the terms of the measurement contribution
					function that are coupled to the connection endpoints */
					if (vs->isEmitterSupernode()) {
						/* If possible, convert 'vt' into an emitter sample */
						if (!vt->cast(m_scene, PathVertex::EEmitterSample) || vt->isDegenerate())
							continue;

						value = radianceWeights[t] *
							vs->eval(m_scene, vsPred, vt, EImportance) *
							vt->eval(m_scene, vtPred, vs, ERadiance);
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
					// 					value *= Path::miWeight(m_scene, m_emitterSubpath, &connectionEdge,
					// 						m_sensorSubpath, s, t, m_sampleDirect, m_lightImage);
					Float weight0 = Path::miWeight(m_scene, m_emitterSubpath, &connectionEdge,
						m_sensorSubpath, s, t, m_sampleDirect, m_lightImage);
					//value *= weight0;
					
					Float miWeight = 0.f;
					MisState sensorState = sensorStates[t - 1];
					MisState emitterState = emitterStates[s - 1];
					if (!sampleDirect){						
						miWeight = miWeightVC(m_scene, vsPred, vs, vt, vtPred,
							s, t,
							emitterState[EVCM], emitterState[EVC],
							sensorState[EVCM], sensorState[EVC],
							misVmWeightFactor, true);
					}else {
						if (t == 1){
							MisState sensorStateTemp;
// 							updateMisHelper(1, m_sensorSubpath, sensorStateTemp, m_scene, 
// 								nLightPaths, misVcWeightFactor, misVmWeightFactor, ERadiance, true);
							PathVertex
								*vtPredTemp = m_sensorSubpath.vertexOrNull(t - 1),
								*vtTemp = m_sensorSubpath.vertex(t);
							miWeight = miWeightVC(m_scene, vsPred, vs, vtTemp, vtPredTemp,
								s, t,
								emitterState[EVCM], emitterState[EVC],
								sensorStateTemp[EVCM], sensorStateTemp[EVC],
								misVmWeightFactor, true);
						}
						else{
							MisState emitterStateTemp;
// 							updateMisHelper(1, m_emitterSubpath, emitterStateTemp, m_scene,
// 								nLightPaths, misVcWeightFactor, misVmWeightFactor, EImportance, true);
							PathVertex
								*vsPredTemp = m_emitterSubpath.vertexOrNull(s - 1),
								*vsTemp = m_emitterSubpath.vertex(s);
							miWeight = miWeightVC(m_scene, vsPredTemp, vsTemp, vt, vtPred,
								s, t,
								emitterStateTemp[EVCM], emitterStateTemp[EVC],
								sensorState[EVCM], sensorState[EVC],
								misVmWeightFactor, true);
						}
					}
					
					value *= miWeight;

					if ((miWeight - weight0) / std::max(miWeight, weight0) > 0.3f){
						float fuck = 1.f;					

						Float weight2 = miWeightVC(m_scene, vsPred, vs, vt, vtPred,
							s, t,
							emitterState[EVCM], emitterState[EVC],
							sensorState[EVCM], sensorState[EVC],
							misVmWeightFactor, true);
						Float weight1 = Path::miWeight(m_scene, m_emitterSubpath, &connectionEdge,
							m_sensorSubpath, s, t, m_sampleDirect, m_lightImage);						

						Float weight3 = miWeightVC(m_scene, vsPred, vs, vt, vtPred,
							s, t,
							emitterState[EVCM], emitterState[EVC],
							sensorState[EVCM], sensorState[EVC],
							misVmWeightFactor, true);
					}

					if (sampleDirect) {
						/* Now undo the previous change */
						if (t == 1)
							m_sensorSubpath.swapEndpoints(vtPred, vtEdge, vt);
						else
							m_emitterSubpath.swapEndpoints(vsPred, vsEdge, vs);
					}

					Spectrum contrib = value;
					if (contrib[0] < 0.f || _isnan(contrib[0]) || !std::isfinite(contrib[0]) ||
						contrib[1] < 0.f || _isnan(contrib[1]) || !std::isfinite(contrib[1]) ||
						contrib[2] < 0.f || _isnan(contrib[2]) || !std::isfinite(contrib[2])){
						Log(EWarn, "Invalid sample value[EPSSMLT-VC]: %f %f %f",
							contrib[0], contrib[1], contrib[2]);
						continue;
					}

					/* Determine the pixel sample position when necessary */
					if (vt->isSensorSample() && !vt->getSamplePosition(vs, samplePos))
						continue;

					if (t < 2) {
						list.append(samplePos, value);
					}
					else {
						BDAssert(m_sensorSubpath.vertexCount() > 2);
						list.accum(0, value);
					}
				}
			}
		}

		std::vector<uint32_t> searchResults;
		std::vector<uint32_t> searchResultsAll;
		std::vector<Point> searchPos;
		std::vector<uint32_t> acceptCnt;
		std::vector<size_t> shootCnt;
		size_t clampThreshold = 1000000; // TODO: make it a parameter
		int shareShootThreshold = 128;

		// massive vertex merging - CAMERA!
		// trace camera subpath to gather photons
		if (useVM){
			PathVertex *vs0 = m_pool.allocVertex();
			PathVertex *vsPred0 = m_pool.allocVertex();
			PathVertex *vsPred20 = m_pool.allocVertex();
			PathVertex *vs = NULL, *vsPred = NULL, *vsPred2 = NULL;
			PathEdge connectionEdge;
			PathVertex *succVertex = m_pool.allocVertex();
			PathEdge *succEdge = m_pool.allocEdge();
			Point2 samplePos(0.0f);			

			int minT = 2; int minS = 3;

			int maxT = (int)m_sensorSubpath.vertexCount() - 1;
			if (m_maxDepth != -1)
				maxT = std::min(maxT, m_maxDepth);
			for (int t = minT; t <= maxT; ++t) {
				PathVertex
					*vtPred2 = m_sensorSubpath.vertex(t - 2),
					*vtPred = m_sensorSubpath.vertex(t - 1),
					*vt = m_sensorSubpath.vertex(t);
				PathEdge
					*predEdge = m_sensorSubpath.edge(t - 2);

				if (!vt->isDegenerate()){
					BDAssert(vt->type == PathVertex::ESurfaceInteraction);

					searchResultsAll.clear();
					m_lightPathTree.search(vt->getPosition(), gatherRadius, searchResultsAll);
					avgGatherPoints.incrementBase();
					avgGatherPoints += searchResultsAll.size();
					maxGatherPoints.recordMaximum(searchResultsAll.size());
					numZeroGatherPoints.incrementBase();
					if (searchResultsAll.size() == 0){
						++numZeroGatherPoints;
						continue;
					}
					
					// select camera direction connection from gather results
					searchResults.clear();
					searchPos.clear();
					Float tBandwidth = 0.f;
					if (vtPred->isSensorSample()){
						tBandwidth = 10000.f;
					}
					else{
						BDAssert(vtPred->isSurfaceInteraction());
						const Intersection &its = vtPred->getIntersection();
						const BSDF *bsdf = its.getBSDF();
						tBandwidth = bsdf->getBandwidth();
					}
					for (int i = 0; i < searchResultsAll.size(); i++){
						LightPathNode node = m_lightPathTree[searchResultsAll[i]];
						size_t vertexIndex = node.data.vertexIndex;
						// decide connection direction
						bool cameraDirConnection = true;
						Float sBandwidth = 0.f;
						vsPred = vsPred0;
						m_lightVerticesExt[vertexIndex - 1].expand(vsPred);						
						if (vsPred->isEmitterSample()){
							PositionSamplingRecord &pRec = vsPred->getPositionSamplingRecord();
							const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
							if (!emitter->needsDirectionSample())
								sBandwidth = 99999.f;
							else
								sBandwidth = 0.f;
						}
						else{
							BDAssert(vsPred->isSurfaceInteraction());
							const Intersection &its = vsPred->getIntersection();
							const BSDF *bsdf = its.getBSDF();
							sBandwidth = bsdf->getBandwidth();
						}
						if (sBandwidth < tBandwidth)
							cameraDirConnection = false;

						if (!cameraDirConnection) continue;

						searchResults.push_back(searchResultsAll[i]);
						LightVertexExt lvertexExt = m_lightVerticesExt[vertexIndex];
						searchPos.push_back(lvertexExt.position);
					}
					acceptCnt.resize(searchResults.size());
					shootCnt.resize(searchResults.size());
					for (int i = 0; i < searchResults.size(); i++){
						acceptCnt[i] = 0;
						shootCnt[i] = 0;
					}

					// evaluate sampling domain pdf normalization
					bool shareShoot = false;
					size_t totalShootShared = 0;
					Float invBrdfIntegralShared = 1.f;
					if (searchPos.size() > shareShootThreshold && tBandwidth < 100 && false){
						Vector4 smplBBox = Vector4(0.f);
						Vector4 smplBBoxDiff = Vector4(0.f);
						Float brdfIntegral = vtPred->gatherAreaPdf(vt->getPosition(), gatherRadius * 2.f, vtPred2, smplBBox, &smplBBoxDiff);
						if (brdfIntegral > 0.f){
							shareShoot = true;
							invBrdfIntegralShared = 1.f / brdfIntegral;
							uint32_t finishCnt = 0;
							Float distSquared = gatherRadius * gatherRadius;
							while (finishCnt < searchResults.size() && totalShootShared < clampThreshold){
								totalShootShared++;

								// restricted sampling evaluation shoots
								if (!vtPred->sampleShoot(m_scene, m_sensorSampler, vtPred2, predEdge, succEdge, succVertex, ERadiance, vt->getPosition(), gatherRadius * 2.f, smplBBox, smplBBoxDiff))
									continue;

								Point pshoot = succVertex->getPosition();
								for (int i = 0; i < searchPos.size(); i++){									
									Float pointDistSquared = (succVertex->getPosition() - searchPos[i]).lengthSquared();
									if (pointDistSquared < distSquared){
										acceptCnt[i]++;
										if (shootCnt[i] == 0)
											finishCnt++;
										shootCnt[i] = totalShootShared;
									}
								}
							}
							avgInvpShoots.incrementBase();
							avgInvpShoots += totalShootShared;
							maxInvpShoots.recordMaximum(totalShootShared);
							numInvpShoots += totalShootShared;
							numClampShoots.incrementBase(searchResults.size());
							if (finishCnt < searchResults.size()){
								numClampShoots += searchResults.size() - finishCnt;
							}
						}
					}					

					MisState sensorState = sensorStates[t - 1];
					for (int k = 0; k < searchResults.size(); k++){
						LightPathNode node = m_lightPathTree[searchResults[k]];
						int s = node.data.depth;
						if (m_maxDepth != -1 && s + t > m_maxDepth + 2 || s < minS) continue;

						if (s == 2 && t == 2) continue;

						size_t vertexIndex = node.data.vertexIndex;
						LightVertex vi = m_lightVertices[vertexIndex];
						LightVertex viPred = m_lightVertices[vertexIndex - 1];
						MisState emitterState = vi.emitterState;

						vs = vs0; vsPred = vsPred0; vsPred2 = vsPred20;
						m_lightVerticesExt[vertexIndex].expand(vs);
						m_lightVerticesExt[vertexIndex - 1].expand(vsPred);
						if (s > 2)
							m_lightVerticesExt[vertexIndex - 2].expand(vsPred2);
						else
							vsPred2->makeEndpoint(m_scene, time, EImportance);

						// decide the direction to do connection
						bool cameraDirConnection = true;
						Float sBandwidth = 0.f, tBandwidth = 0.f;
						if (vtPred->isSensorSample()){
							tBandwidth = 10000.f;
						}
						else{
							BDAssert(vtPred->isSurfaceInteraction());
							const Intersection &its = vtPred->getIntersection();
							const BSDF *bsdf = its.getBSDF();
							tBandwidth = bsdf->getBandwidth();
						}
						if (vsPred->isEmitterSample()){
							PositionSamplingRecord &pRec = vsPred->getPositionSamplingRecord();
							const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
							if (!emitter->needsDirectionSample())
								sBandwidth = 99999.f;
							else
								sBandwidth = 0.f;
						}
						else{
							BDAssert(vsPred->isSurfaceInteraction());
							const Intersection &its = vsPred->getIntersection();
							const BSDF *bsdf = its.getBSDF();
							sBandwidth = bsdf->getBandwidth();
						}
						if (sBandwidth < tBandwidth)
							cameraDirConnection = false;

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
							int interactions = m_maxDepth - s - t + 1;
							if (contrib.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vs, vtPred, NULL, interactions))
								continue;
							contrib *= connectionEdge.evalCached(vs, vtPred, PathEdge::EGeneralizedGeometricTerm);
						}
						else{
							contrib = radianceWeights[t] * viPred.importanceWeight * invLightPaths;
							contrib *= vt->eval(m_scene, vtPred, vsPred, ERadiance) *	vsPred->eval(m_scene, vsPred2, vt, EImportance);
							int interactions = m_maxDepth - s - t + 1;
							if (contrib.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vt, vsPred, NULL, interactions))
								continue;
							contrib *= connectionEdge.evalCached(vt, vsPred, PathEdge::EGeneralizedGeometricTerm);
						}
						
						Float invp = 0.f;						
						if (acceptCnt[k] > 0){
							invp = (Float)(shootCnt[k]) / (Float)(acceptCnt[k])* invBrdfIntegralShared;
						}else{
							Vector4 smplBBox = Vector4(0.f);
							Vector4 smplBBoxDiff = Vector4(0.f);
							Float brdfIntegral;
							if (cameraDirConnection)
								brdfIntegral = vtPred->gatherAreaPdf(vs->getPosition(), gatherRadius, vtPred2, smplBBox, &smplBBoxDiff);
							else
								brdfIntegral = vsPred->gatherAreaPdf(vt->getPosition(), gatherRadius, vsPred2, smplBBox, &smplBBoxDiff);

							if (brdfIntegral == 0.f) continue;
							Float invBrdfIntegral = 1.f / brdfIntegral;
							size_t totalShoot = 0, acceptedShoot = 0, targetShoot = 1;
							Float distSquared = gatherRadius * gatherRadius;
							while (totalShoot < clampThreshold){
								totalShoot++;

								// restricted sampling evaluation shoots
								Float pointDistSquared;
								if (cameraDirConnection){
									if (!vtPred->sampleShoot(m_scene, m_lightPathSampler, vtPred2, predEdge, succEdge, succVertex, ERadiance, vs->getPosition(), gatherRadius, smplBBox, smplBBoxDiff))
										continue;
									pointDistSquared = (succVertex->getPosition() - vs->getPosition()).lengthSquared();
								}
								else{
									if (!vsPred->sampleShoot(m_scene, m_lightPathSampler, vsPred2, predEdge, succEdge, succVertex, EImportance, vt->getPosition(), gatherRadius, smplBBox, smplBBoxDiff))
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
							avgInvpShoots.incrementBase();
							avgInvpShoots += totalShoot;
							maxInvpShoots.recordMaximum(totalShoot);
							numInvpShoots += totalShoot;
							numClampShoots.incrementBase();
							if (totalShoot >= clampThreshold)
								++numClampShoots;
							invp = (acceptedShoot > 0) ? (Float)(totalShoot) / (Float)(acceptedShoot)* invBrdfIntegral : 0;
						}						
						contrib *= invp;

						// accumulate to image
						if (contrib.isZero()) continue;

						// MIS weighting						
						Float miWeight = miWeightVM(m_scene, s, t, vsPred2, vsPred, vs, vt, vtPred, vtPred2,
							emitterState, sensorState, misVcWeightFactor, cameraDirConnection);

						contrib *= miWeight;

						if (contrib[0] < 0.f || _isnan(contrib[0]) || !std::isfinite(contrib[0]) ||
							contrib[1] < 0.f || _isnan(contrib[1]) || !std::isfinite(contrib[1]) ||
							contrib[2] < 0.f || _isnan(contrib[2]) || !std::isfinite(contrib[2])){
							Log(EWarn, "Invalid sample value[EPSSMLT-UPM]: %f %f %f, invp = %f, miWeight = %f",
								contrib[0], contrib[1], contrib[2], invp, miWeight);
								//contrib[0], contrib[1], contrib[2], invpOrig, miWeight);								
							continue;
						}

						if (t == 2) {
							list.append(samplePos, contrib);
						}
						else {
							list.accum(0, contrib);
						}
					}
				}
			}
			m_pool.release(vs0);
			m_pool.release(vsPred0);
			m_pool.release(vsPred20);
			m_pool.release(succVertex);
			m_pool.release(succEdge);
		}

		// massive vertex merging - LIGHT!
		// trace light subpath to gather importons
		//int maxImporton = 100;
		if (useVM){
			PathVertex *vt0 = m_pool.allocVertex();
			PathVertex *vtPred0 = m_pool.allocVertex();
			PathVertex *vtPred20 = m_pool.allocVertex();
			PathVertex *vt = NULL, *vtPred = NULL, *vtPred2 = NULL;
			PathEdge connectionEdge;
			PathVertex *succVertex = m_pool.allocVertex();
			PathEdge *succEdge = m_pool.allocEdge();
			Point2 samplePos(0.0f);
			std::vector<uint32_t> searchResults;

			int minT = 2; int minS = 3;

			int maxS = (int)m_emitterSubpath.vertexCount() - 1;
			if (m_maxDepth != -1)
				maxS = std::min(maxS, m_maxDepth);
			for (int s = minS; s <= maxS; ++s) {
				PathVertex
					*vsPred2 = m_emitterSubpath.vertex(s - 2),
					*vsPred = m_emitterSubpath.vertex(s - 1),
					*vs = m_emitterSubpath.vertex(s);
				PathEdge
					*predEdge = m_emitterSubpath.edge(s - 2);

				if (!vs->isDegenerate()){
					BDAssert(vs->type == PathVertex::ESurfaceInteraction);

					searchResultsAll.clear();
					m_cameraPathTree.search(vs->getPosition(), gatherRadius, searchResultsAll);					

					// select camera direction connection from gather results
					searchResults.clear();
					searchPos.clear();
					Float sBandwidth = 0.f;
					if (vsPred->isEmitterSample()){
						PositionSamplingRecord &pRec = vsPred->getPositionSamplingRecord();
						const Emitter *emitter = static_cast<const Emitter *>(pRec.object);
						if (!emitter->needsDirectionSample())
							sBandwidth = 99999.f;
						else
							sBandwidth = 0.f;
					}
					else{
						BDAssert(vsPred->isSurfaceInteraction());
						const Intersection &its = vsPred->getIntersection();
						const BSDF *bsdf = its.getBSDF();
						sBandwidth = bsdf->getBandwidth();
					}					
					for (int i = 0; i < searchResultsAll.size(); i++){
						LightPathNode node = m_cameraPathTree[searchResultsAll[i]];
						size_t vertexIndex = node.data.vertexIndex;						
						int t = node.data.depth;
						if (s == 2 && t == 2) continue;

// 						LightVertexExt lvertexExt = m_cameraVerticesExt[vertexIndex];
// 						if (dot(vsPred->getPosition() - lvertexExt.position, lvertexExt.geoFrameN) < 0.f) continue;

						// decide connection direction
						bool cameraDirConnection = true;
						vtPred = vtPred0;
						m_cameraVerticesExt[vertexIndex - 1].expand(vtPred);
						Float tBandwidth = 0.f;
						if (vtPred->isSensorSample()){
							tBandwidth = 10000.f;
						}
						else{
							BDAssert(vtPred->isSurfaceInteraction());
							const Intersection &its = vtPred->getIntersection();
							const BSDF *bsdf = its.getBSDF();
							tBandwidth = bsdf->getBandwidth();
						}
						if (sBandwidth < tBandwidth)
							cameraDirConnection = false;

						if (cameraDirConnection) continue;
						LightVertexExt lvertexExt = m_cameraVerticesExt[vertexIndex];
						searchResults.push_back(searchResultsAll[i]);						
						searchPos.push_back(lvertexExt.position);
					}
					acceptCnt.resize(searchResults.size());
					shootCnt.resize(searchResults.size());
					for (int i = 0; i < searchResults.size(); i++){
						acceptCnt[i] = 0;
						shootCnt[i] = 0;
					}
					avgGatherPoints.incrementBase();
					avgGatherPoints += searchResultsAll.size();
					maxGatherPoints.recordMaximum(searchResultsAll.size());
					numZeroGatherPoints.incrementBase();
					if (searchResultsAll.size() == 0){
						++numZeroGatherPoints;
						continue;
					}

					// evaluate sampling domain pdf normalization
					float confidence = 0.95f;
					Float expectShoot = log(1.f - pow(confidence, 1.f / float(searchPos.size()))) / log(0.75f);
					bool shareShoot = false;
					size_t totalShootShared = 0;
					Float invBrdfIntegralShared = 1.f;
					if (searchPos.size() > expectShoot){
						Vector4 smplBBox = Vector4(0.f);
						Vector4 smplBBoxDiff = Vector4(0.f);
						Float brdfIntegral = vsPred->gatherAreaPdf(vs->getPosition(), gatherRadius * 2.f, vsPred2, smplBBox, &smplBBoxDiff);
						if (brdfIntegral > 0.f){
							shareShoot = true;
							invBrdfIntegralShared = 1.f / brdfIntegral;
							uint32_t finishCnt = 0;
							Float distSquared = gatherRadius * gatherRadius;
							//while (finishCnt < searchResults.size() && totalShootShared < clampThreshold){
							while (finishCnt < searchResults.size() && totalShootShared < expectShoot){
								totalShootShared++;

								// restricted sampling evaluation shoots
								if (!vsPred->sampleShoot(m_scene, m_lightPathSampler, vsPred2, predEdge, succEdge, succVertex, EImportance, vs->getPosition(), gatherRadius * 2.f, smplBBox, smplBBoxDiff))
									continue;

								Point pshoot = succVertex->getPosition();
								for (int i = 0; i < searchPos.size(); i++){
									Float pointDistSquared = (succVertex->getPosition() - searchPos[i]).lengthSquared();
									if (pointDistSquared < distSquared){
										acceptCnt[i]++;
										if (shootCnt[i] == 0)
											finishCnt++;
										shootCnt[i] = totalShootShared;
									}
								}
							}
							avgInvpShoots.incrementBase();
							avgInvpShoots += totalShootShared;
							maxInvpShoots.recordMaximum(totalShootShared);
							numInvpShoots += totalShootShared;
							numClampShoots.incrementBase(searchResults.size());
							if (finishCnt < searchResults.size()){
								numClampShoots += searchResults.size() - finishCnt;
							}
						}
					}

					//Float keepProb = (float)maxImporton / (float)searchResults.size();
					MisState emitterState = emitterStates[s - 1];
					for (int k = 0; k < searchResults.size(); k++){
						LightPathNode node = m_cameraPathTree[searchResults[k]];
						int t = node.data.depth;
						if (m_maxDepth != -1 && s + t > m_maxDepth + 2 || t < minT) continue;

						if (s == 2 && t == 2) continue;

						//if (m_lightPathSampler->next1D() > keepProb) continue;
						
						size_t vertexIndex = node.data.vertexIndex;
						LightVertex vi = m_cameraVertices[vertexIndex];
						LightVertex viPred = m_cameraVertices[vertexIndex - 1];
						MisState sensorState = vi.emitterState;

						vt = vt0; vtPred = vtPred0; vtPred2 = vtPred20;
						m_cameraVerticesExt[vertexIndex].expand(vt);
						m_cameraVerticesExt[vertexIndex - 1].expand(vtPred);
						if (t > 2)
							m_cameraVerticesExt[vertexIndex - 2].expand(vtPred2);
						else
							vtPred2->makeEndpoint(m_scene, time, ERadiance);

						// all light direction due to pre selection
						bool cameraDirConnection = false;

						// get image position from importon
						samplePos.x = vi.wo.x;
						samplePos.y = vi.wo.y;

						// evaluate contribution
						Spectrum contrib;
						contrib = importanceWeights[s - 1] * vi.importanceWeight * invLightPaths;
						contrib *= vt->eval(m_scene, vtPred, vsPred, ERadiance) *	vsPred->eval(m_scene, vsPred2, vt, EImportance);
						int interactions = m_maxDepth - s - t + 1;
						if (contrib.isZero() || !connectionEdge.pathConnectAndCollapse(m_scene, NULL, vt, vsPred, NULL, interactions))
							continue;
						contrib *= connectionEdge.evalCached(vt, vsPred, PathEdge::EGeneralizedGeometricTerm);

						//contrib /= keepProb;

						// unbiased evaluate connection probability 1/p 
						Float invp = 0.f;
						if (acceptCnt[k] > 0){
							invp = (Float)(shootCnt[k]) / (Float)(acceptCnt[k])* invBrdfIntegralShared;
						}
						else{
							Vector4 smplBBox = Vector4(0.f);
							Vector4 smplBBoxDiff = Vector4(0.f);
							Float brdfIntegral = vsPred->gatherAreaPdf(vt->getPosition(), gatherRadius, vsPred2, smplBBox, &smplBBoxDiff);
							if (brdfIntegral == 0.f) continue;
							Float invBrdfIntegral = 1.f / brdfIntegral;
							size_t totalShoot = 0, acceptedShoot = 0, targetShoot = 1;
							Float distSquared = gatherRadius * gatherRadius;
							while (totalShoot < clampThreshold){
								totalShoot++;

								// restricted sampling evaluation shoots
								Float pointDistSquared;
								if (!vsPred->sampleShoot(m_scene, m_lightPathSampler, vsPred2, predEdge, succEdge, succVertex, EImportance, vt->getPosition(), gatherRadius, smplBBox, smplBBoxDiff))
									continue;
								pointDistSquared = (succVertex->getPosition() - vt->getPosition()).lengthSquared();

								if (pointDistSquared < distSquared){
									acceptedShoot++;
									if (acceptedShoot == targetShoot){
										break;
									}
								}
							}
							avgInvpShoots.incrementBase();
							avgInvpShoots += totalShoot;
							maxInvpShoots.recordMaximum(totalShoot);
							numInvpShoots += totalShoot;
							numClampShoots.incrementBase();
							if (totalShoot >= clampThreshold)
								++numClampShoots;
							invp = (acceptedShoot > 0) ? (Float)(totalShoot) / (Float)(acceptedShoot)* invBrdfIntegral : 0;
						}
						contrib *= invp;

						// accumulate to image
						if (contrib.isZero()) continue;

						// MIS weighting						
						Float miWeight = miWeightVM(m_scene, s, t, vsPred2, vsPred, vs, vt, vtPred, vtPred2,
							emitterState, sensorState, misVcWeightFactor, cameraDirConnection);

						contrib *= miWeight;

						if (contrib[0] < 0.f || _isnan(contrib[0]) || !std::isfinite(contrib[0]) ||
							contrib[1] < 0.f || _isnan(contrib[1]) || !std::isfinite(contrib[1]) ||
							contrib[2] < 0.f || _isnan(contrib[2]) || !std::isfinite(contrib[2])){
							Log(EWarn, "Invalid sample value[EPSSMLT-UPM]: %f %f %f, invp = %f, miWeight = %f",
								contrib[0], contrib[1], contrib[2], invp, miWeight);
							//contrib[0], contrib[1], contrib[2], invpOrig, miWeight);								
							continue;
						}

						if (t == 2 || true) {
							list.append(samplePos, contrib);
						}
						else {
							list.accum(0, contrib);
						}
					}
				}
			}
			m_pool.release(vt0);
			m_pool.release(vtPred0);
			m_pool.release(vtPred20);
			m_pool.release(succVertex);
			m_pool.release(succEdge);
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

		list.append(samplePos, value);
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
MTS_IMPLEMENT_CLASS(UPMWorkResult, false, WorkResult)
MTS_NAMESPACE_END
