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
		return new ImageBlock(Bitmap::ESpectrum,
			m_film->getCropSize(), m_film->getReconstructionFilter());
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
			true, m_config.initialRadius, m_sampler);
	}

	void process(const WorkUnit *workUnit, WorkResult *workResult, const bool &stop) {		
		ImageBlock *result = static_cast<ImageBlock *>(workResult);
		ImageBlock *midres = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize(), m_film->getReconstructionFilter());
		midres->clear();
		const SeedWorkUnit *wu = static_cast<const SeedWorkUnit *>(workUnit);
		const int workID = wu->getID();
		SplatList *splats = new SplatList();		

		HilbertCurve2D<int> hilbertCurve;
		TVector2<int> filmSize(m_film->getCropSize());
		hilbertCurve.initialize(filmSize);
		uint64_t iteration = workID;
		size_t actualSampleCount = 0;

		float radius = m_config.initialRadius;
		ref<Timer> timer = new Timer();
		for (size_t j = 0; j < m_config.sampleCount || (wu->getTimeout() > 0 && (int)timer->getMilliseconds() < wu->getTimeout()); j++) {
			if (m_config.initialRadius > 0.0f){
				Float reduceFactor = 1.f / std::pow(Float(iteration + 1), (Float)(0.5 * (1 - 0.75/*radiusAlpha*/)));
				radius = std::max(reduceFactor * m_config.initialRadius, (Float)1e-7);
				iteration += 8;
			}
			m_pathSampler->gatherLightPaths(m_config.useVC, m_config.useVM, radius, hilbertCurve.getPointCount(), midres);

			for (size_t i = 0; i < hilbertCurve.getPointCount(); ++i) {
				if (stop) break;

				Point2i offset = Point2i(hilbertCurve[i]);
				m_sampler->generate(offset);
				m_pathSampler->sampleSplatsVCM(m_config.useVC, m_config.useVM, radius, offset, i, *splats);

				for (size_t k = 0; k < splats->size(); ++k) {
					Spectrum value = splats->getValue(k);
					midres->put(splats->getPosition(k), &value[0]);
				}
			}
			actualSampleCount++;
		}

		Log(EInfo, "Run %d iterations", actualSampleCount);

		const Spectrum *pmidres = (Spectrum *)midres->getBitmap()->getData();
		Spectrum *presult = (Spectrum *)result->getBitmap()->getData();
		size_t pixelCount = midres->getBitmap()->getPixelCount();
		for (int i = 0; i < pixelCount; i++){
			presult[i] = pmidres[i] / float(actualSampleCount);;
		}
		midres->clear();
		//delete midres; // TODO: better handle the memory of midres
		
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
	size_t pixelCount = m_accum->getBitmap()->getPixelCount();
	const Spectrum *accum = (Spectrum *) m_accum->getBitmap()->getData();
	Spectrum *target = (Spectrum *) m_developBuffer->getData();

	for (size_t i=0; i<pixelCount; ++i) {
		target[i] = accum[i] / float(m_config.workUnits);
	}
	m_film->setBitmap(m_developBuffer);
	m_refreshTimer->reset();

	m_queue->signalRefresh(m_job);
}

void VCMProcess::processResult(const WorkResult *wr, bool cancelled) {
	LockGuard lock(m_resultMutex);
	const ImageBlock *result = static_cast<const ImageBlock *>(wr);
	m_accum->put(result);
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
		m_accum = new ImageBlock(Bitmap::ESpectrum, m_film->getCropSize());
		m_accum->clear();
		m_developBuffer = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, m_film->getCropSize());
	}
}

MTS_IMPLEMENT_CLASS_S(VCMRenderer, false, WorkProcessor)
MTS_IMPLEMENT_CLASS(VCMProcess, false, ParallelProcess)
MTS_IMPLEMENT_CLASS(SeedWorkUnit, false, WorkUnit)

MTS_NAMESPACE_END
