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

#include "../LibImportanceTypes.h"
#include "Particle.h"
#include "Hit.h"
#include "PointKdTreeShared.h"
#include "..\caching\CacheStats.h"
#include "Stack.h"
#include "Utils.h"
#include "particlefilter.h"

#pragma warning(push)
#pragma warning(disable:4127)

namespace Importance {

    struct EmInfo;

    /************************************************************************/
    /*    DataFilter - filters out particles from opposite hemisphere       */
    /************************************************************************/

    /** 
     *  Filters out particles from opposite hemisphere and also turns their direction. 
     *  Filtering can be switched off by setting isDistOverSphere to true.
     *  Input of constructor is a pointer to particle map (vector). 
     *  Then prepareData gets indices of N-nearest particles in the map.
     */
    template<class TDataVector> 
    class DataFilter {
    public:
        DataFilter( const IStack<Importance::Particle> * photonMap, const Config & config ) :
          m_photonMap( photonMap ), m_cfg( config ), m_info( NULL ) {
              IMPORTANCE_ASSERT( photonMap != NULL );
          }

          inline void prepareData (
              const KdQueryResult         * searchResults,
              const Importance::Hit       & hit,
              int                           nPhotonsFound,
              CacheStats                  & cstats,
              TDataVector                 & dataset ) {

                  int count = 0;
                  Float sqrSearchRadius = 0.f;
                  for( int i = 0; i < nPhotonsFound; ++i ) {			    
                      IMPORTANCE_ASSERT( searchResults[ i ].index >= 0 && searchResults[ i ].index < m_photonMap->size() );
                      Importance::Particle photon   = (*m_photonMap)[ (size_t) searchResults[ i ].index ];
                      photon.incidentDir            = -photon.incidentDir;	

                      if ( ParticleFilterPredicate()( photon.incidentDir, photon.normal, hit.normal, m_cfg.fitting.isFitOverSphere ) ) {
                              ++ count;
                              sqrSearchRadius = std::max( sqrSearchRadius, searchResults[ i ].distSqr );
                              //IMPORTANCE_ASSERT( photon.weight > 0 );
                              dataset.add( photon, hit.normal, photon.weight );
                              if ( m_info ) {
                                  m_info->photonsQueryRes.push( searchResults[ i ].index );
                                  Importance::Particle p = photon;
                                  p.incidentDir = -p.incidentDir;
                                  m_info->particles.push( p );
                              }
                      }
                  }
                  cstats.particleSearchCount    = count;
                  cstats.sqrSearchRadius        = sqrSearchRadius;
          }

          /* If clamp factor is greater than zero it clamps weights of particles by median value times weightClampFactor.
             In either case it computes statistics about particle values used for training distribution dist.
          */
          static inline void clampAndStats( TDataVector & dataset, EmInfo * info, CacheStats & cacheStats, Float weightClampFactor ) {
              float photonDist      = 0.f,
                    photonWeightSum = 0.f,              
                    max             = 0, 
                    min             = std::numeric_limits<float>::max(),
                    maxDist         = 0.f,
                    threshold       = 0.f;
              
              // compute clamping threshold
              bool isClamp = weightClampFactor > 0.f && dataset.size() > 2;
              if ( isClamp ) {
                std::nth_element( dataset.begin(), dataset.begin() + dataset.size() / 2, dataset.end(), dataset.getWeightComparator() );
                threshold = dataset[dataset.size()/2].getSingleValue() * weightClampFactor;
              }


              for( int i = 0; i < dataset.size(); i++ ) {                  
                  Float weight          = dataset[ i ].getSingleValue();
                  Float dist            = dataset[ i ].distance; 
                  //IMPORTANCE_ASSERT( weight > 0.0f );
                  IMPORTANCE_ASSERT( isReal( weight ) );

                  //if clamping is allowed and necessary then clamp
                  if ( isClamp && threshold < weight ) {
                      if ( info != NULL ) {
                        ++info->clampedWeights;
                      }
                      dataset[ i ].setValue( threshold );
                      weight = dataset[ i ].getSingleValue();
                  }

                  IMPORTANCE_ASSERT( isReal( dataset[ i ].distance ) );

                  photonWeightSum += weight;
                  photonDist += weight * dataset[ i ].distance;                  

                  if ( weight < min ) { min = weight; }
                  if ( weight > max ) { max = weight; }
                  if ( dist > maxDist ) { maxDist = dist; }
              }

              cacheStats.photonWeightSum += photonWeightSum;
              cacheStats.photonDist += photonDist;
              cacheStats.photonCount += (int) dataset.size();
              cacheStats.maxDistance = std::max( std::sqrt( maxDist ), cacheStats.sqrSearchRadius );
              cacheStats.avgParticleDistance = std::max(cacheStats.photonDist, 1e-5f)/std::max(1e-5f, cacheStats.photonWeightSum);
              IMPORTANCE_ASSERT(isReal(cacheStats.avgParticleDistance) && cacheStats.avgParticleDistance >= 0);

              if ( info != NULL ) {
                info->minWeight = std::min( min, info->minWeight );
                info->maxWeight = std::max( max, info->maxWeight );
                info->avgWeight = cacheStats.photonWeightSum / cacheStats.photonCount;
              }
          }

          inline void setInfo( EmInfo * info ) {
              m_info = info;
          }
    private:
        void operator=( const DataFilter<TDataVector> & ) {}
    private:
        EmInfo * m_info;
        const IStack<Importance::Particle> * m_photonMap;
        const Config & m_cfg;
    };


    /************************************************************************/
    /* Abstract distribution factory                                        */
    /************************************************************************/

    struct Hit;
    struct Config;
    class DefaultDistributionModel;


    /************************************************************************/
    /* BasicDistributionFactory */
    /************************************************************************/

    class BasicDistributionFactory {
    public:     
        virtual ~BasicDistributionFactory() { }

        virtual void init ( const IStack<Importance::Particle>& particles, const Importance::Config * config ) {
            m_photonMap     = &particles;
            m_config        = config;             
        }

        virtual bool getInstance( 
            const Hit& hit, 
            const KdQueryResult* particles,
            int nParticlesFound, 
            DefaultDistributionModel & outDist ) {
                throw std::runtime_error( "Not implemented" );
        }

        virtual bool getInstanceProgressive( 
            const Hit& /*hit*/, 
            const KdQueryResult* /*particles*/,
            int /*nParticlesFound*/, 
            DefaultDistributionModel & /*outDist*/ ) {
                throw std::runtime_error( "Not implemented" );
        }

        virtual const IStack<Importance::Particle> * getPhotonMap() const {
            return m_photonMap;
        }
        virtual const Config * getConfig() const {
            return m_config;
        }

        IMPORTANCE_INLINE void computeGradient( const Importance::Hit & hit, KdQueryResult * particles, const int nParticlesFound, CacheStats & cs, float & outGradient  )
        {
            // UPDATE gradient
            Vector3 xDir, yDir;
            coordinateSystem(hit.normal, xDir, yDir);
            IMPORTANCE_ASSERT(xDir.isNormalized() && yDir.isNormalized());
            
            for(int i = 0; i < nParticlesFound; ++i) {
                const Particle& particle = getPhotonMap()->get(particles[i].index);
                if(particle.weight > 0.f && dot(particle.incidentDir, hit.normal) < 0.f) {
                    //const float distance = (particle.position-hit.position).size();
                    cs.irradiance += particle.weight;
                    const Vector3 normDir = (particle.position-hit.position).getNormalized();
                    cs.dX += particle.weight*dot(normDir, xDir);
                    cs.dY += particle.weight*dot(normDir, yDir);
                }
            }

            const Float area = PI * cs.maxDistance * cs.maxDistance;
            const Float d_irradiance = cs.irradiance / area;
            const Float d_dX = cs.dX / ( area * cs.maxDistance );
            const Float d_dY = cs.dY / ( area * cs.maxDistance );

            outGradient = Vector3(d_dX/std::max(1e-6f, d_irradiance), d_dY/std::max(1e-6f, d_irradiance), 0).size();
            //IMPORTANCE_ASSERT(outGradient.isValid());
        }

        virtual std::string toString() const { return "toString() method is not implemented yet."; }
    public:
        const IStack<Importance::Particle> * m_photonMap;
        const Importance::Config * m_config;
    };

}

#pragma warning(pop)