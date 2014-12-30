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

#include <mitsuba/render/libImpVizUtils.h>
#include <mitsuba/mitsuba.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/render/libImpUtils.h>
#include <mitsuba/render/weightwindow.h>

MTS_NAMESPACE_BEGIN 

/** Holds and trains samplers for guiding both paths from light sources and from the camera. */
class GuidingSamplers : public ConfigurableObject {
public:

	GuidingSamplers::GuidingSamplers() : ConfigurableObject(Properties()) {
		m_timer = new Timer();
		m_trainingTimer = new Timer();
		m_failureListMutex = new Mutex();
	}

	GuidingSamplers(const Properties &props)
		  : ConfigurableObject(props), m_cfg(GuidingConfig(props)),
            m_qmcSamplerID_photons( -1 ), m_qmcSamplerID_importons( -1 ), 
            m_photonTracingState( 0 ), m_importonTracingState( 0 ),
            m_radianceSampler( NULL ), m_importanceSampler( NULL ), m_enviroSampler( NULL ),
            m_scene( NULL ), m_canceled( false ) {
        m_timer					= new Timer();
        m_trainingTimer			= new Timer();
        m_failureListMutex		= new Mutex();
    }

	ref<GuidingSamplers> clone(){
		ref<GuidingSamplers> sampler = new GuidingSamplers();
		sampler->m_cfg = m_cfg;
		sampler->m_qmcSamplerID_photons = m_qmcSamplerID_photons;
		sampler->m_qmcSamplerID_importons = m_qmcSamplerID_importons;
		sampler->m_photonTracingState = m_photonTracingState;
		sampler->m_importonTracingState = m_importonTracingState;
		sampler->m_radianceSampler = m_radianceSampler;
		sampler->m_importanceSampler = m_importanceSampler;
		sampler->m_enviroSampler = m_enviroSampler;
		sampler->m_scene = m_scene;
		sampler->m_canceled = false;
		
		sampler->m_radianceStats.reset();
		sampler->m_importanceStats.reset();
		sampler->m_enviroSamplerStats.reset();
		sampler->m_failureList.clear();

		sampler->m_ww = m_ww;
		return sampler.get();
	}

    void trainingPhase( const RenderJob *job, int sceneResID, int sensorResID ) {
        if ( m_cfg.m_mitsuba.useGuidedSampling) {
            if ( !m_cfg.m_mitsuba.usePingPong ) {
                SLog( EInfo, "Ping-pong is switched off - setting number of training passes to 1." );
                m_cfg.m_mitsuba.nPasses = 1;
            }
        }
        else {
            SLog( EInfo, "Path guiding is switched off - setting number of training passes to 0." );
            m_cfg.m_mitsuba.nPasses = 0;
        }

        // Starts training phase
        m_trainingTimer->reset();
        for ( int iPass = 0; iPass < m_cfg.m_mitsuba.nPasses && !m_canceled; ++iPass ) {
            if ( m_cfg.m_mitsuba.usePingPong || ( m_cfg.m_mitsuba.useEnvironmentSampler && iPass == 0 ) ) {	
                // Evaluates as true for all passes if ping pong is used or for first pass if only env. sampler is used
                /* IMPORTONS TRACING */
                ref<BgParticlesVec> bgParticles;
                SLog( EInfo, "Training (pass %i/%i): tracing importons", iPass+1, m_cfg.m_mitsuba.nPasses);
                ref<PhotonMap> importonsMap = particleTracingPass( job, sceneResID, sensorResID, m_radianceSampler, GuidedEmitParticleProcess::EFromSensor, &bgParticles );
                SLog( EInfo, "Training: importance cache update" );
                if ( !m_importanceSampler) {
                    m_importanceSampler = Importance::SamplerFactory::createSampler( m_cfg.m_importance );
                    SAssert(m_importanceSampler);
                    m_importanceSampler->init( PhotonsIterator( importonsMap ),
                        m_cfg.m_importance,
                        getImportanceCamera(), &m_importanceStats );
                }
                else {
                    m_importanceSampler->refreshSamples( PhotonsIterator( importonsMap ) );
                }
                /* EM SAMPLER CREATE & UPDATE */
                if ( m_cfg.m_mitsuba.useEnvironmentSampler ) {
                    // Find environment sampler
                    const EnviroMapInterface * enviroMap = dynamic_cast<const EnviroMapInterface *>(m_scene->getEnvironmentEmitter());
                    if ( enviroMap == NULL ) {
                        SLog( EWarn, "The environment light source was not found." );
                        m_cfg.m_mitsuba.useEnvironmentSampler = false;
                    } else {
                        if ( !m_enviroSampler ) {
                            SLog( EInfo, "Training: creating environment sampler");													
                            m_enviroSampler = Importance::SamplerFactory::createEnviroSampler( m_cfg.m_importance );
                            SAssert(m_enviroSampler);
                            enviroMap->setEnviroSampler(m_enviroSampler);
                            bool isOk = m_enviroSampler->init( EnviroImportonsIterator( bgParticles ),
                                m_cfg.m_importance,
                                &m_enviroSamplerStats,
                                enviroMap );
                            if ( !isOk ) {
                                SLog( EWarn, "Initialization of environment sampler has failed. \n %s", m_enviroSampler->toString().c_str() );
                                delete m_enviroSampler;
                                m_enviroSampler = NULL;
                                m_cfg.m_mitsuba.useEnvironmentSampler = false;
                                enviroMap->setEnviroSampler(NULL);
                            } 
                        }
                        else {
                            SLog( EInfo, "Training (pass %i/%i): training environment sampler", iPass+1, m_cfg.m_mitsuba.nPasses );
                            m_enviroSampler->refreshSamples( EnviroImportonsIterator( bgParticles ), enviroMap );
                        }
                    }
//#ifdef LIBIMP_STATS
//                    std::ostringstream fname;
//                    fname << "importanceDirectionsPdf_" << std::setw( 2 ) << std::setfill( '0' ) << iPass << ".exr";
//                    writeImportanceEnviroMap( fname.str(), m_enviroSampler->getDebugMaps()->m_pdf );
//
//                    fname = std::ostringstream();
//                    fname << "importanceDirectionsProduct_" << std::setw( 2 ) << 
//                        std::setfill( '0' ) << iPass << ".exr";
//                    writeImportanceEnviroMap( fname.str(), *m_enviroSampler->getDirectionsBmp() );
//#endif
                }                
                SLog( EInfo, "%s", cacheStatsToString( m_importanceStats, "Importance" ).c_str() );
            }
            /* PHOTONS TRACING */
            SLog( EInfo, "Training (pass %i/%i): tracing photons", iPass+1, m_cfg.m_mitsuba.nPasses );
            ref<PhotonMap> photonMap = particleTracingPass( job, sceneResID, sensorResID, m_importanceSampler, GuidedEmitParticleProcess::EFromEmitters, NULL );
            SLog( EInfo, "Training: radiance cache update" );
            if ( !m_radianceSampler) {
                m_radianceSampler = Importance::SamplerFactory::createSampler( m_cfg.m_importance );
                SAssert(m_radianceSampler);
                m_radianceSampler->init( PhotonsIterator( photonMap ),
                    m_cfg.m_importance,
                    getImportanceCamera(), &m_radianceStats );
            }
            else{
                m_radianceSampler->refreshSamples( PhotonsIterator( photonMap ) );
            }
            SLog( EInfo, "%s", cacheStatsToString( m_radianceStats, "Radiance" ).c_str() );
        } // end of the main training loop


        if ( m_canceled ) {
            SLog( EWarn, "The training phase was stopped by the user." );
        }

        if ( m_cfg.m_mitsuba.useGuidedSampling ) {
            SLog( EInfo, "Training phase took %s", timeString( m_trainingTimer->getMilliseconds() / 1000.f ).c_str() );
        }                
    }



    bool preprocess( const Scene * scene ) {
        m_scene    = scene;
        m_canceled = false;        

        if ( !m_cfg.m_mitsuba.useGuidedSampling ) {
            return true;
        }

        Importance::initialize();

        m_cfg.computeBBox( scene );

        m_importonTracingState = m_photonTracingState = 0;

        m_radianceStats.reset();
        m_importanceStats.reset();
        m_enviroSamplerStats.reset();
        m_failureList.clear();        

        m_ww = WeightWindow( m_cfg.m_mitsuba.m_ww.lowerBound * 1e-6f, 
                             m_cfg.m_mitsuba.m_ww.size,
                             Spectrum( m_scene->emitterPdfSum() ) /* Approx. total flux in the scene. */,
                             Spectrum( 1.f ) /* Sensor's "importance power" */ );

        return true;
    }

	void done(){
		if (m_qmcSamplerID_photons > -1) {
			m_qmcSamplerID_photons = -1;
		}

		if (m_qmcSamplerID_importons > -1) {
			m_qmcSamplerID_importons = -1;
		}

		m_radianceSampler = NULL;
		m_importanceSampler = NULL;
		if (m_enviroSampler) {
			m_enviroSampler = NULL;
		}
		m_failureList.clear();
	}

    void postprocess() {
        if ( !m_cfg.m_mitsuba.useGuidedSampling ) {
            return;
        }             

        statistics_and_visual( m_scene, m_cfg, m_radianceSampler, NULL, &m_radianceStats, 
            &m_enviroSamplerStats, false, &m_failureList );

        if ( m_cfg.m_mitsuba.usePingPong ) {
            statistics_and_visual( m_scene, m_cfg, m_importanceSampler, m_enviroSampler, 
                &m_importanceStats, &m_enviroSamplerStats, true, &m_failureList );
        }

        if ( m_qmcSamplerID_photons > -1 ) {
            Scheduler::getInstance()->unregisterResource( m_qmcSamplerID_photons );
            m_qmcSamplerID_photons = -1;
        }

        if ( m_qmcSamplerID_importons > -1 ) {
            Scheduler::getInstance()->unregisterResource( m_qmcSamplerID_importons );
            m_qmcSamplerID_importons = -1;
        }

        Importance::SamplerFactory::disposeSampler( m_radianceSampler );
        m_radianceSampler = NULL;
        Importance::SamplerFactory::disposeSampler( m_importanceSampler );
        m_importanceSampler = NULL;
        if ( m_enviroSampler ) {
            delete m_enviroSampler; 
            m_enviroSampler = NULL;
            const EnviroMapInterface * enviroMap = dynamic_cast<const EnviroMapInterface *>(m_scene->getEnvironmentEmitter());
            enviroMap->setEnviroSampler( NULL );
        }

        m_failureList.clear();

        Importance::close();
    }

    /** User has canceled the process. */
    void cancel() {
        m_canceled = true;
    }

    const WeightWindow & getWeightWindow() const {
        return m_ww;
    }

    WeightWindow & getWeightWindow() {
        return m_ww;
    }

    void distributionConstructionFailed( const Intersection& its ) const {
#ifdef LIBIMP_STATS
        m_failureListMutex->lock();
        m_failureList.push_back( itsToHit( its ));
        m_failureListMutex->unlock();
#endif
    }


    std::string toString() const {
        std::ostringstream oss;
        if ( m_cfg.m_importance.cache.useCache ) {
            oss << cacheStatsToString( m_importanceStats, "Importance" );
        }        

        if ( m_cfg.m_importance.cache.useCache && m_cfg.m_mitsuba.usePingPong ) {
            oss << cacheStatsToString( m_radianceStats, "Radiance" );
        }

        return oss.str();
    }


    ~GuidingSamplers() {
		SLog(EInfo, "Destruct %d", (void*)this);
//         SAssert( m_enviroSampler == NULL );
//         SAssert( m_importanceSampler == NULL );
//         SAssert( m_radianceSampler == NULL );
//         SAssert( m_qmcSamplerID_photons == -1 && m_qmcSamplerID_importons == -1 );
    }

    /************************************************************************/
    /* Getters                                                              */
    /************************************************************************/
    Importance::Sampler * getRadianceSampler() const { return m_radianceSampler; }
    Importance::Sampler * getImportanceSampler() const { return m_importanceSampler; }
    Importance::IEnviroSampler * getEnviroSampler() const { return m_enviroSampler; }
    Importance::Sampler * getRadianceSampler() { return m_radianceSampler; }
    Importance::Sampler * getImportanceSampler() { return m_importanceSampler; }
    Importance::IEnviroSampler * getEnviroSampler() { return m_enviroSampler; }
    const GuidingConfig & getConfig() const { return m_cfg; }


protected:

	/** Traces particles either from camera or from light sources, importance can either be null or 
	should store radiance sampler for tracing importons or importance sampler for tracing photons */
	ref<PhotonMap> particleTracingPass( const RenderJob * job, int sceneResID,             
		int sensorResID, Importance::Sampler * importance,
		GuidedEmitParticleProcess::ETraceDirection traceDirection,
		/* out */ ref<BgParticlesVec> * bgImportons ) {

			m_timer->reset();

			// Get tracing direction specific settings
			int & samplerID = traceDirection == GuidedEmitParticleProcess::EFromSensor ? m_qmcSamplerID_importons : m_qmcSamplerID_photons;
			size_t & tracingState = traceDirection == GuidedEmitParticleProcess::EFromSensor ? m_importonTracingState : m_photonTracingState;
			const size_t & nParticles = traceDirection == GuidedEmitParticleProcess::EFromSensor ? m_cfg.m_mitsuba.nImportons : m_cfg.m_mitsuba.nPhotons;
            std::string particleName( traceDirection == GuidedEmitParticleProcess::EFromSensor ? 
                "importon" : "photon" );
			if ( samplerID == -1 ) { /** Creates sampler if it does not exists */
				SLog( EInfo, "Creating independent samplers for %s tracing, one sampler per core...", particleName.c_str() );
				ref<Scheduler> sched = Scheduler::getInstance();
				ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(Sampler), Properties(/*"halton"*/"independent")));
				/* Create a sampler instance for every core */
				size_t unCoreCount = sched->getCoreCount();
				std::vector<SerializableObject *> samplers(sched->getCoreCount());
				for (size_t i=0; i < unCoreCount ; ++i) {
					ref<Sampler> clonedSampler = sampler->clone();
					clonedSampler->incRef();
					samplers[i] = clonedSampler.get();
				}
				samplerID = sched->registerMultiResource(samplers);
				tracingState = 0;
			}

            if ( traceDirection == GuidedEmitParticleProcess::EFromSensor ) {
                getWeightWindow().pathTracing();
            } else {
                getWeightWindow().lightTracing();
            }

			/* Generate the photon map */
			ref<Scheduler> sched = Scheduler::getInstance();
			ref<GuidedEmitParticleProcess> proc = new GuidedEmitParticleProcess(
    			GuidedEmitParticleProcess::EAllIndirectSurfaceParticles,                
                //GuidedGatherPhotonProcess::EAllSurfacePhotons,
				nParticles,
				/*m_granularity*/0, 
				m_cfg.m_mitsuba.maxDepth == -1 ? -1 : m_cfg.m_mitsuba.maxDepth - 1, 
				0,
                getWeightWindow(),
				true,
				true, 
				NULL,
				importance,				
				traceDirection,
                m_cfg,
				tracingState );

			proc->bindResource("scene", sceneResID);
			proc->bindResource("sensor", sensorResID);
			proc->bindResource("sampler", samplerID);

			sched->schedule(proc);
			sched->wait(proc);
			tracingState = proc->getNumGenerated();

			if (proc->getReturnStatus() != ParallelProcess::ESuccess) {
				return NULL;
			}

			SLog( EInfo, "Global photon map full. Shot " SIZE_T_FMT " particles, excess particles due to parallelism: " 
				SIZE_T_FMT " Stored particles: " SIZE_T_FMT, proc->getShotParticles(), proc->getExcessPhotons(), proc->getPhotonMap()->size() );
			SLog( EInfo, "Tracing took %s", timeString( m_timer->getMilliseconds() / 1000.f ).c_str() );

			if ( bgImportons ) {
				*bgImportons = proc->getBgParticles();
				SLog( EInfo, "Captured " SIZE_T_FMT " background particles.", 
					(*bgImportons)->m_data.size() );
			}

            SLog( EInfo, "Weight window: %s", getConfig().m_mitsuba.useWeightWindow ? "ON" : "OFF" );
            SLog( EInfo, "Weight window: \n%s", getWeightWindow().toString().c_str() );
            SLog( EInfo, "Particle (%s) statistics:\n%s", particleName.c_str(), proc->getParticleStats().toString().c_str() );
            
			return proc->getPhotonMap();
	}


	/** Returns a camera compatible with ImportanceLib */
	inline const Importance::Camera * getImportanceCamera() {
		return new MitsubaImportanceCamera( m_scene->getSensor() );
	}  

private:
    /** Illumination importance function used for guiding indirect rays */
    Importance::Sampler * m_radianceSampler;
    /** Illumination importance function used for guiding photons */
    Importance::Sampler * m_importanceSampler;
    /** Importance library sampler of environment map */
    Importance::IEnviroSampler * m_enviroSampler;


    /** List with hits were distributions could not be created (Debugging and visualization) */
    mutable std::vector<Importance::Hit> m_failureList;
    /** Mutex - used for changing a failure list */
    mutable ref<Mutex> m_failureListMutex;

    /** Sampler ID for photons and importons tracing */
    int m_qmcSamplerID_photons, m_qmcSamplerID_importons;
    /** Photons and importons tracing state */
    size_t m_photonTracingState, m_importonTracingState;

    /** Statistics */
    Importance::Stats m_radianceStats;
    Importance::Stats m_importanceStats;
    Importance::Stats m_enviroSamplerStats;

    /** Has user canceled rendering? */
    bool m_canceled;

    /** Configuration */
    GuidingConfig m_cfg;

    /** Weight window */
    mitsuba::WeightWindow m_ww;
    
    const Scene * m_scene;

    /** Measurement of training phase */
    ref<Timer> m_trainingTimer;
    /** Timer for measurement of single particle tracing steps. */
    ref<Timer> m_timer;
};

MTS_NAMESPACE_END
