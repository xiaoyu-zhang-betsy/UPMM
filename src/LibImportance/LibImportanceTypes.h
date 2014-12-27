/*
    This file is part of LibImportance library that provides a technique for guiding
    transport paths towards the important places in the scene. This is a direct implementation
    of the method described in the paper "On-line Learning of Parametric Mixture 
    Models for Light Transport Simulation", ACM Trans. Graph. (SIGGRAPH 2014) 33, 4 (2014).
   
    Copyright (c) 2014 by Jiri Vorba, Ondrej Karlik, Martin Sik.

    LibImportance library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    LibImportance library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once

#include <vector>

#include "shared/Vector2.h"
#include "shared/Vector3.h"
#include "shared/Utils.h"
#include "shared/StaticArray.h"
#include "shared/BoundingBox.h"
#include "shared/PointKdTreeShared.h"
#include "shared/simplelogger.h"
#include "shared/Particle.h"
#include "shared/Hit.h"
#include <iostream>
#include "shared/Stack.h"

#pragma warning(push)
#pragma warning(disable:4100)
#pragma warning(disable:4127)

namespace Importance {    


    /************************************************************************/
    /*  Ancestor of all result buffers */
    /************************************************************************/

    struct IResultBuffer {};


    /************************************************************************/
    /*  Distribution - user interface of importance lib                     */
    /************************************************************************/

    struct EmInfo;

    class Distribution {
    public:

        /// The user interface

        virtual void sampleDirections(const Vector2* random, Vector3 * output, float* pdfs, const int count) = 0;

        virtual void pdfs(const Vector3 * directions, float* pdfs, const int count) const = 0;

        IMPORTANCE_INLINE Float pdf( const Vector3 & direction ) {
            Float res;
            pdfs( &direction, &res, 1 );
            return res;
        }

        virtual void release() = 0;

        virtual std::string toString() const {
            return "Distribution: toString() is not implemented";
        }

        virtual ~Distribution() { }
        

        //////////////////////////////////////////////////////////////////////////
        /// Developement stuff, DANGEROUS stuff *x*

        virtual const EmInfo * getInfo() const {
            return NULL;
        }

        virtual const void * getModel() const {
            return NULL;
        }
    };

    /************************************************************************/
    /*  Information about fitting process                                   */
    /************************************************************************/

    struct EmInfo {
        EmInfo() {
            clear();
        };

        void clear() {
            nPhotonsUsed        = 0;
            constructionFail    = false;
            lkhFunction.clear();
            nPhotonsIgnored     = 0;
            photonsQueryRes.clear();
            nIteration          = 0;
            photonmap           = NULL;
            timeEM              = 0;
            timeInit            = 0;
            nPhotonsFound       = 0;
            nnsearchRadius      = 0.f;
			clampedWeights		= 0;
			maxWeight			= -std::numeric_limits< float >::max();
			minWeight			= std::numeric_limits< float >::max();
			avgWeight			= 0;
            particles.clear();
        }
        void release() {
            lkhFunction.dealloc();
            photonsQueryRes.dealloc();
            particles.dealloc();
        }

        size_t nPhotonsUsed;
        bool constructionFail;
        IStack< Float > lkhFunction;        
        size_t nPhotonsIgnored;
        IStack< int > photonsQueryRes;
        const IStack< Importance::Particle > * photonmap;
        IStack< Particle > particles;
        size_t nIteration;
        unsigned int timeEM;
        unsigned int timeInit;
        size_t nPhotonsFound;
        Float nnsearchRadius;
		unsigned int clampedWeights;
		float maxWeight;
		float minWeight;
		float avgWeight;

        void merge( const EmInfo & info ) {            
            nPhotonsUsed = info.nPhotonsUsed;
            constructionFail = info.constructionFail;
lkhFunction = info.lkhFunction;            
            nPhotonsIgnored = info.nPhotonsIgnored;
            photonsQueryRes = info.photonsQueryRes;
            photonmap = info.photonmap;
            nIteration = info.nIteration;

            /// add particles
            for(auto it = info.particles.cbegin(); it != info.particles.cend(); ++it) {
                particles.push(*it);
            }

            timeEM = info.timeEM;
            timeInit = info.timeInit;
            nPhotonsFound = info.nPhotonsFound;
            nnsearchRadius =info.nnsearchRadius;
        }
    };


    /************************************************************************/
    /* DistributionType enum                                                */
    /************************************************************************/

    enum EDistributionType {        
        EVmfMixture,
        EPharr,
        EHey,			
        EGaussianMixture,
        EJensen
    };

    /************************************************************************/
    /* RegularizationType enum                                                */
    /************************************************************************/

    enum ERegularizationType {
        ENone,
        EMAP,
        ESatoIshy
    };

    /* Type of environment importance distributions (if it is used) */
    enum EEnviroType {
        EDynamic,
        EStatic,
        EStaticDirections,
        EStaticPositions
    };

    /************************************************************************/
    /* Camera                                                               */
    /************************************************************************/

    class Camera {
    public:
        virtual Float pixelToWorldRatio(const Hit& hit) const = 0;
    };


    /************************************************************************/
    /* Config                                                               */
    /************************************************************************/

    struct Config {        

        struct Cache {
            /** invGeometryError from Ward's irradiance caching scheme. We set this to one, so it has no effect at the moment. */
            float invGeometryError;
    
            /** Min validity radius */
			float minRadius;

            /** Max validity radius */
            float maxRadius;

            /**
              Max environment records search radius for both particles and cache records
              - in solid degrees              
              */
            float maxEnviroRadius;

            /** Kullback-Leibler divergence threshold */
            float klDivThres;

            /** Multiply of particle search radius that is an upper bound for validity radius */
            float photonClampingMax;

           /** Multiply of particle search radius that is a lower bound for validity radius */
            float photonClampingMin;

            /** Neighbour clamping of cache records on/off */
            bool neighbourClamping;

            /** Caching - on/off */
            bool useCache;
                        
            /** 
               An average number of observed particles per mixture component that must be used 
               for distribution training to consider the trained distribution to be reliable. 
               If the distribution is not reliable then the validity radius is not determined 
               based on KL divergence criterion.

               This is used in static stepwise EM.
            */
            unsigned int undersampledCoef;

            /** 
               Same as undersampledCoef.            

               This is used in on-line stepwise EM.
            */
            unsigned int undersampledCoefOnline;



            Cache() {
				useCache                    = true;
                minRadius                   = 0.1f;
                maxRadius                   = 1.f;
                /* Environment search radius is given in multiplies of 1 degree */
                maxEnviroRadius             = 1.f;                
                invGeometryError            = 1.f;
                klDivThres       = 5.0f;
                neighbourClamping           = true;
                photonClampingMax          =  1.f;
                photonClampingMin           = 0.5f;
                undersampledCoef            = 5;
                undersampledCoefOnline      = undersampledCoef;     
            }
        } cache;

        // if true, camera pixel to world is used, else scene bbox size times globalSizeMult
        bool useCameraPixel2World;
        /* This influence the maximum size of search radii of both particles and cache records. 
           These are set relatively to the scene size. This global parameter is set in multiples 
           of scene size. */
        float globalSizeMult;        

        struct Particles {
            int knn;
            Particles() {
                knn = 250;
            }
        } particles;

        struct VMFM {
            enum InitType {
                EKMeansPP,
                ERandom,
                ESingular,
                EPhotonUniform,
                EPhotonWeighed,
                ERandRestart
            } init;            
            size_t nPhotonsPerComponent;
            Float maxKappa;
            size_t nRandomRestarts;             

            /// developement - screen resolution for PerPixelCacheSampler
            Vector2i res;

            VMFM() {
                nPhotonsPerComponent = 0;
                /// -1.f means no user restriction
                maxKappa    = 0.5e6f;

                init        = ERandom;       
                                
                res         = Vector2i( 512, 512 );

                /// implicitly there are no random restarts
                nRandomRestarts         = 1;
            }            
        } vmfm;

		struct Gaussian
		{			
			Float lambda;
			size_t rsemIter;
			Gaussian() : lambda(0.1f), rsemIter(10){}
		} gaussian;

        struct Pharr {
            Pharr() {
                cosAngleDir = std::cos( 10.f / 180.f * IMP_PI );
            }
            Float cosAngleDir;
        } pharr;      

        struct Jensen {
            Jensen() {
                w = h = 16;
            }
            // grid resolution
            int w, h;
        } jensen;

        ELogLevel logLevel;

        Vector3f sceneBboxMin;
        Vector3f sceneBboxMax;

        bool perfectSampling;

        //developement
        bool isWeightOverride;

        // What kind of distribution to use
        EDistributionType distType;
        
        struct MisCorrection {
            bool use;
            float vmfSamplingFraction;

            MisCorrection() {
                use = false;
                vmfSamplingFraction = NAN;
            }
        } misCorrection;

        struct Fitting {
            Fitting() {
                 /// how many components use for the distribution based on number of used photons in knn
                 nMaxComponents                 = 8;
                 nComponentsEnviroDirections    = 16;
                 /// start without static EM on the first batch
                 isPureOnline       = false;
                 /// batch size for online EM (how often do we update the model)
                 batchSize          = 10;
                 /// forgetting parameter for online EM
                 alpha              = 0.7f;
                 /// parameter for online EM (number of photons after which the algorithm updates the model for the first time )
                 delayedUpdate      = 1; 
                 // gather statistics about fitting process
                 isUseStats         = true;
                 isFitOverSphere    = false;
                 regularizationType = EMAP;
                 /// maximal number of EM iteration 
                 nMaxEMIter          = 250;
				 weightClampFactor	 = -1.f;
                 enviroWeightClampFactor = -1.f;
                 // number of random restarts
                 nRandomRestarts     = 0; 
				 // number of converged records divided by total records which is neccessary to stop progressive updating
				 convergedRecordsRatio = 0.8f;
                 enviroType            = EDynamic;
            }

            ERegularizationType regularizationType;
            /* Works for every factory (Hey, Pharr, Vmfm, ...) */
            bool isFitOverSphere;
            bool isPureOnline;            
            int batchSize;
            unsigned int delayedUpdate;
            Float alpha;            
            bool isUseStats;
			Float weightClampFactor;
            Float enviroWeightClampFactor;
            size_t nMaxEMIter;
            unsigned int nRandomRestarts;
			float convergedRecordsRatio;
            EEnviroType enviroType;
            size_t nMaxComponents;
            size_t nComponentsEnviroDirections;

        } fitting;

        Config() {
            perfectSampling       = false;
            logLevel              = EDebug;
            distType              = EGaussianMixture;
            useCameraPixel2World  = false;
            globalSizeMult        = 1/100.f;
            isWeightOverride      = false;
        }

        IMPORTANCE_INLINE float pixel2WorldRatio(const Hit& hit, const Camera* camera) const {
            return useCameraPixel2World ? camera->pixelToWorldRatio(hit) : (sceneBboxMax-sceneBboxMin).size()*globalSizeMult;
        }

        
    };

    /************************************************************************/
    /*  Input iterator                                                     */
    /************************************************************************/


    class InputIterator {
    public:
        virtual ~InputIterator() { }

        virtual Vector3f position() const = 0;

        virtual Vector3f incidentDirection() const = 0;

        virtual float weight() const = 0;

        /// \brief The distance the Photon travelled before hitting the surface
        virtual float distance() const = 0;

        /// Normal in the photon's hitpoint
        virtual Vector3f normal() const = 0;

        virtual void next() = 0;

        virtual bool isValid() const = 0;
    };

    template<class TStatsInt64>
    struct StatsBase {

        struct {
            TStatsInt64 attempted;
            TStatsInt64 rejected;
        } queries;
        
        struct {
            TStatsInt64 records;
            TStatsInt64 recordsFound;          
            TStatsInt64 newRecords;
            TStatsInt64 newUndersampledRecords;
        } cache;

        void reset() {
            memset(this, 0, sizeof(*this));
        }

        void nextIteration() {
            cache.newRecords              = 0;
            cache.newUndersampledRecords  = 0;
        }
    };

    typedef StatsBase<StatsInt64> Stats;
};

#pragma warning(pop)