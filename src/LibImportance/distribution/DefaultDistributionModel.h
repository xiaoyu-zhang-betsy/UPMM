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

#include "../shared/basicfactory.h"

namespace Importance {

    /************************************************************************/
    /* Common ancestor of various distributions used for fitting            */
    /************************************************************************/        

    class DefaultDistributionModel {
    public:
        DefaultDistributionModel() {
            m_undersampled = false;
        }

        IMPORTANCE_INLINE static bool getDistribution( BasicDistributionFactory* factory, 
            const Hit& hit, 
            KdQueryResult * particles, 
            const int nParticlesFound, 
            DefaultDistributionModel& outDist ) {            
                const bool res = factory->getInstance(hit, particles, nParticlesFound, outDist );
                IMPORTANCE_ASSERT( res );            
                return res;
        }

        IMPORTANCE_INLINE static bool getDitributionUpdated( BasicDistributionFactory * factory, 
            const Importance::Hit & hit, 
            KdQueryResult * particles, 
            const int nParticlesFound, 
            DefaultDistributionModel & outDist ) { 			
                const bool res = factory->getInstanceProgressive( hit, particles, nParticlesFound, outDist );
                return res;                
        }

        virtual void release() {}
        virtual const EmInfo * getInfo() const { return NULL; }
        virtual std::string toString() const { return "toString() method is not implemented"; }         

        //////////////////////////////////////////////////////////////////////////
        // Caching stuff - this must implement every distribution which can be placed into cache.
        //////////////////////////////////////////////////////////////////////////
        static IMPORTANCE_INLINE Float validityRadius(const DefaultDistributionModel& dist, const Config& config, const float pixel2world = 0.f ) {
            //IMPORTANCE_ASSERT( false );
            //return -1.f;
            return 1.f;
        }

        template<class TRecord>
        static IMPORTANCE_INLINE Float weight(const Hit& hit, const TRecord& record, const Config& config) {
            //IMPORTANCE_ASSERT( false );
            //return -1.f;
            return 1.f;
        }        

        virtual ~DefaultDistributionModel() { }

        IMPORTANCE_INLINE bool isUndersampled() const { return m_undersampled; }
        IMPORTANCE_INLINE void setUndersampled( bool val ) { m_undersampled = val; }

        IMPORTANCE_INLINE float getPhotonMaxDistance() const {
            return m_cacheStats.maxDistance;
        }		

    public:
        /// Caching statistics 
        CacheStats  m_cacheStats;
    private:
        /// information for caching - says whether we can trust the distribution
        bool m_undersampled;
    };
}