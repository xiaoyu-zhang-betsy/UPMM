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

#include "../LibImportanceTypes.h"

#include "../pharr/pharrfactory.h"
#include "../hey/heyfactory.h"
#include "../gaussian/gaussian_implsse.h"
#include "../jensen/jensenfactory.h"
#include "../enviro/EnviroSampler.h"
#include "../shared/buffers.h"

#pragma warning (disable: 4127)

namespace Importance {

    /************************************************************************/
    /* Typedefs                                                             */
    /************************************************************************/

    typedef SimpleSampler<Jensen>               JensenSimpleSampler;
    typedef SimpleSampler<Pharr>                PharrSimpleSampler;
    typedef SimpleSampler<Hey>                  HeySimpleSampler;
	typedef SimpleSampler<GaussianMixtureSSEHemisphere>	GaussianMixtureSimpleSampler;


   // typedef CachedSampler<VMFM>                         VMFMCachedSampler;
	typedef CachedSampler<GaussianMixtureSSEHemisphere> GaussianMixtureCachedSampler;
    typedef CachedSampler<Pharr>                        PharrCachedSampler;
    typedef CachedSampler<Hey>                          HeyCachedSampler;
    typedef CachedSampler<Jensen>                       JensenCachedSampler;

    /************************************************************************/
    /* SamplerFactory                                                       */
    /************************************************************************/

    class SamplerFactory {
    public:
        static IEnviroSampler * createEnviroSampler( const Config & cfg ) {
            IEnviroSampler * sampler = NULL;

            switch ( cfg.distType ) {
            case EPharr:
                sampler = new EnviroSampler<Pharr>(); 
                break;
            case EHey:
                sampler = new EnviroSampler<Hey>();
                break;
            case EGaussianMixture:
                sampler = new EnviroSampler<GaussianMixtureSSEHemisphere>();
                break;
            case EJensen:
                sampler = new EnviroSampler<Jensen>();
                break;
            }

            sampler->setDistributionFactory( createDistributionFactory( cfg.distType ) );
            return sampler;
        }

        static Sampler * createSampler( const Config & cfg ) {
            Sampler * sampler = NULL;

            switch ( cfg.distType ) {
                case EPharr:
                    if ( cfg.cache.useCache ) {
                        sampler = new Importance::PharrCachedSampler();
                    } else {                    
                        sampler = new Importance::PharrSimpleSampler(); 
                    }
                    break;
                case EHey:
                    if ( cfg.cache.useCache ) {
                        sampler = new Importance::HeyCachedSampler();
                    } else {
                        sampler = new Importance::HeySimpleSampler();
                    }
                    break;
                case EGaussianMixture:
                    if ( cfg.cache.useCache ) {
                        sampler = new Importance::GaussianMixtureCachedSampler();
                    } else {
                        sampler = new Importance::GaussianMixtureSimpleSampler();
                    }
                    break;
                case EJensen:
                    if ( cfg.cache.useCache ) {
                        sampler = new Importance::JensenCachedSampler();
                    } else {
                        sampler = new Importance::JensenSimpleSampler();
                    }
                    break;
            }

            sampler->setDistributionFactory( createDistributionFactory( cfg.distType ) );
            return sampler;
        }

        static void disposeSampler( Sampler * sampler ) {
            if ( sampler ) { delete sampler; }
        }

   private:
        /** Creates the right distribution factory based on a distribution type */
        static BasicDistributionFactory * createDistributionFactory( EDistributionType distType ) {
            switch( distType ) {
                case EPharr:
                    return new PharrFactory();
                case EHey:
                    return new HeyFactory();
                case EGaussianMixture:
                    return new GaussianMixtureFactorySSEHemisphere();
                case EJensen:
                    return new JensenFactory();
            }

            IMPORTANCE_ASSERT( false );
            return NULL;
        }
    };
}