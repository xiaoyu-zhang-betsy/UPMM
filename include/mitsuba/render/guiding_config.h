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
#include "LibImportance/LibImportance.h"

MTS_NAMESPACE_BEGIN

struct GuidingConfig {   
	GuidingConfig() {}
    GuidingConfig( const Properties & props ) {   
        // Load Importance Library properties

		// Cache
		m_importance.cache.minRadius				= props.getFloat("minRadius",m_importance.cache.minRadius);
		m_importance.cache.maxRadius				= props.getFloat("maxRadius",m_importance.cache.maxRadius);
		m_importance.cache.maxEnviroRadius			= props.getFloat("maxEnviroRadius",m_importance.cache.maxEnviroRadius);
		m_importance.cache.klDivThres	= props.getFloat("klDivThres",m_importance.cache.klDivThres);
		m_importance.cache.photonClampingMin		= props.getFloat("photonClampingMin",m_importance.cache.photonClampingMin);
		m_importance.cache.photonClampingMax		= props.getFloat("photonClampingMax",m_importance.cache.photonClampingMax);
		m_importance.cache.neighbourClamping		= props.getBoolean("neighbourClamping",m_importance.cache.neighbourClamping);
		m_importance.cache.useCache					= props.getBoolean("useCache",m_importance.cache.useCache);

		// KNN
		m_importance.particles.knn					= props.getInteger("knnParticles",m_importance.particles.knn);

		// Camera global multiplication
		m_importance.globalSizeMult					= props.getFloat( "globalSizeMult", m_importance.globalSizeMult);

		// Jensen method settings
		m_importance.jensen.w		                = props.getInteger( "jensenRes", m_importance.jensen.w );
		m_importance.jensen.h						= m_importance.jensen.w;

		// Distribution to be fit
		std::string iifName							= props.getString("IIF", "gaussian");

		if ( iifName == "pharr" ) {
            m_importance.distType = Importance::EPharr;
        } else if ( iifName == "hey" ) {
            m_importance.distType = Importance::EHey;
        } else if ( iifName == "gaussian" ) {
            m_importance.distType = Importance::EGaussianMixture;
        } else if ( iifName == "jensen" ) {
            m_importance.distType = Importance::EJensen;
        } else {                
            SLog( EError, "Unknown importance function \"%s\".", iifName.c_str() );
        }

		// Load other properties (related only to mitsuba)

		m_mitsuba.nPasses					= props.getInteger( "passes", (int)m_mitsuba.nPasses );
		m_mitsuba.nPhotons					= props.getInteger( "nPhotons", (int)m_mitsuba.nPhotons );
		m_mitsuba.nImportons				= props.getInteger( "nImportons", (int)m_mitsuba.nImportons );
		m_mitsuba.showVisualization			= props.getBoolean( "showVisualization", m_mitsuba.showVisualization );
		m_mitsuba.usePingPong				= props.getBoolean( "usePingPong", m_mitsuba.usePingPong );
		m_mitsuba.useEnvironmentSampler		= props.getBoolean( "useEnvSampler", m_mitsuba.useEnvironmentSampler );
		m_mitsuba.useGuidedSampling			= props.getBoolean( "useGuidedSampling", m_mitsuba.useGuidedSampling );
		m_mitsuba.bsdfSamplingProbability	= props.getFloat( "bsdfSamplingProbability", m_mitsuba.bsdfSamplingProbability );
        m_mitsuba.maxDepth                  = props.getInteger( "maxDepth", m_mitsuba.maxDepth );
        m_mitsuba.useWeightWindow           = props.getBoolean( "useWeightWindow", m_mitsuba.useWeightWindow );
    }

	inline void computeBBox( const Scene * scene ) {
		/// bounding box
		const Point & p1 = scene->getAABB().min;
		const Point & p2 = scene->getAABB().max;
		m_importance.sceneBboxMin.x = p1.x;
		m_importance.sceneBboxMin.y = p1.y;
		m_importance.sceneBboxMin.z = p1.z;
		m_importance.sceneBboxMax.x = p2.x;
		m_importance.sceneBboxMax.y = p2.y;
		m_importance.sceneBboxMax.z = p2.z;
	}


    std::string samplerTypeStr(const Importance::EDistributionType & distType, bool suffix = false) const {
        std::string res = "";
        switch( distType ) {
        case Importance::EPharr:
            res = "EPharr";
            break;
        case Importance::EHey:
            res = "EHey";
            break;
        case Importance::EGaussianMixture:
            res = "EGaussianMixture";
            break;
        case Importance::EJensen:
            res = "EJensen";
            break;
        default:
            res = "unknown";
            SLog( EError, "Unknown smapler type" );                
        }
		if (!suffix) {
			return res;
        } else {
			return res + (m_importance.cache.useCache ? " - CACHED" : " - WITHOUT CACHE");
        }
    }

	/** Properties to string */
	std::string toString() const {
		std::ostringstream oss;
		// Distribution to be fit
		oss	<< "  minRadius \t\t= " << m_importance.cache.minRadius << "," << std::endl
			<< "  maxRadius \t\t= " << m_importance.cache.maxRadius << "," << std::endl
			<< "  maxEnviroRadius \t\t= " << m_importance.cache.maxEnviroRadius << "," << std::endl
			<< "  klDivThres \t\t= " << m_importance.cache.klDivThres << "," << std::endl
			<< "  photonClampingMin \t\t= " << m_importance.cache.photonClampingMin << "," << std::endl
			<< "  photonClampingMax \t\t= " << m_importance.cache.photonClampingMax << "," << std::endl
			<< "  neighbourClamping \t\t= " << m_importance.cache.neighbourClamping << "," << std::endl
			<< "  useCache \t\t\t= " << m_importance.cache.useCache << "," << std::endl
			<< "  knnParticles \t\t= " << m_importance.particles.knn << "," << std::endl
			<< "  globalSizeMult \t\t= " << m_importance.globalSizeMult << "," << std::endl
			<< "  jensenRes \t\t= " << m_importance.jensen.w << "," << std::endl
			<< "  IFF \t\t\t= " << samplerTypeStr(m_importance.distType,true) << "," << std::endl
			<< "  passes \t\t\t= " << m_mitsuba.nPasses << "," << std::endl
			<< "  nPhotons \t\t\t= " << m_mitsuba.nPhotons << "," << std::endl
			<< "  nImportons \t\t= " << m_mitsuba.nImportons << "," << std::endl
			<< "  showVisualization \t\t= " << m_mitsuba.showVisualization << "," << std::endl
			<< "  usePingPong \t\t= " << m_mitsuba.usePingPong << "," << std::endl
			<< "  useEnvSampler \t\t= " << m_mitsuba.useEnvironmentSampler << "," << std::endl
			<< "  useGuidedSampling \t\t= " << m_mitsuba.useGuidedSampling << "," << std::endl
			<< "  bsdfSamplingProbability \t= " << m_mitsuba.bsdfSamplingProbability << std::endl;
		return oss.str();
	}

	Importance::Config m_importance;

	struct MitsubaConfig {
		/** Number of emitted photons*/
		size_t nPhotons;
		/** Number of emitted importons*/
		size_t nImportons;
		/** Show visualization after rendering?*/
		bool showVisualization;
		/** Use ping-pong?*/
		bool usePingPong;
		/** Number of training passes*/
		size_t nPasses;
		/** Use guided environment sampler?*/
		bool useEnvironmentSampler;
		/** Bsdf sampling probability */
		Float bsdfSamplingProbability;
		/** Guided sampling? */
		bool useGuidedSampling;
        /** Maximum lenght of a light-path (number of segments) */
        int maxDepth;
        /** Weight window on/off */
        bool useWeightWindow;

        struct WeightWindowConfig {
            /** Lower bound of weight window in multiplies of 1e-6.
                If the weight window is used the particles with weight under the lowerBound
                are killed with probability particleWeight/lowerBound.
                This ensures that particle weights never drop under the lowerBound.
            */
            Float lowerBound;
            /** Size of the weight window (i.e. the weight window upper bound = lowerBound * 1e-6 + size).
                Particles with the weight above the upper bound are split. */
            Float size;

            WeightWindowConfig() {
                lowerBound  = 1.f;
                size        = 2.f;
            }
        } m_ww;

		MitsubaConfig() {
			nPhotons                = 100000;
			nImportons              = 100000;
			showVisualization       = false;
			usePingPong             = true;
			nPasses                 = 10;
			useEnvironmentSampler   = false;
			useGuidedSampling       = true;
			bsdfSamplingProbability = 0.5f;
            maxDepth                = -1;
            useWeightWindow         = true;
		}
	};

	MitsubaConfig m_mitsuba;

};

MTS_NAMESPACE_END