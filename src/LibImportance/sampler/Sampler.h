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
#include "../shared/basicfactory.h"

namespace Importance {    

    /************************************************************************/
    /*  Sampler                                                             */
    /************************************************************************/

    class VizAPI;

    class Sampler {
    protected:

        virtual Distribution* getDistributionImpl(const Hit& hit, IResultBuffer& buffer) = 0;
        virtual void initImpl(InputIterator& samples, const Config& config, const Camera* camera, Stats* stats) = 0;

    public:

        BasicDistributionFactory * distFactory;

        Stats* stats;

    public:
        Sampler() : distFactory( NULL ) {}

        virtual ~Sampler() { 
            if ( distFactory != NULL ) { delete distFactory; }            
        }   

        void setDistributionFactory( BasicDistributionFactory * factory ) { distFactory = factory; }

        IMPORTANCE_INLINE void init(InputIterator& samples, const Config& config, const Camera* camera, Stats* stats) {
            this->stats = stats;
            //memset(stats, 0, sizeof(*stats));
            initImpl(samples, config, camera, stats);
            SimpleLogger::getInstance()->setLogLevel( config.logLevel );
        }

        virtual bool isCacheUsed() const = 0;

        virtual VizAPI * getVizApi() = 0;

#ifdef LIBIMP_GATHER_PARTICLES
        virtual void saveGatherPoints( const std::vector<Hit> & hits ) = 0;
        virtual void loadGatherPoints() = 0;
#endif

        /// \brief Returns pointer to created distribution (int the buffer memory), or NULL, 
        IMPORTANCE_INLINE Distribution* getDistribution(const Hit& hit, IResultBuffer& buffer) {
            IMPORTANCE_INC_STATS(stats->queries.attempted);
            Distribution* distr = getDistributionImpl(hit, buffer);
            if(distr == NULL) {
                IMPORTANCE_INC_STATS(stats->queries.rejected);
            }
            return distr;
        }

        /// Implement in your sampler to utilize progressivity
        virtual void refreshSamples( InputIterator & samples, bool refineCache = true ) {
            throw std::runtime_error( "refreshSamples() method is not implemented" );
        }


        virtual std::string toString() const { return "Method toString() is not implemented."; };

        ///developement 
        virtual bool setCreateNewRecords( bool ) { return false; }
    };
}