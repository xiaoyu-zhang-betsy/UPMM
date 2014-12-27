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

#include "State.h"
#include "..\shared\basicfactory.h"

namespace Importance {
    
    template<typename TDistributionModel>
    class SamplersFacade{    
    public:
        SamplersFacade( IEnviroSampler * enviroSampler, SimpleSampler<TDistributionModel> * sampler ) :
          m_enviroSampler( enviroSampler ), m_sampler( sampler ) {}

        IMPORTANCE_INLINE SimpleSampler<TDistributionModel> * getSampler() { return m_sampler; }

        IMPORTANCE_INLINE const SimpleSampler<TDistributionModel> * getSampler() const { return m_sampler; }

        IMPORTANCE_INLINE IEnviroSampler * getEnviroSampler() { return m_enviroSampler; }

        IMPORTANCE_INLINE const IEnviroSampler * getEnviroSampler() const { return m_enviroSampler; }

        int nnquery( const Hit & hit, ImStaticArray<KdQueryResult,MAX_KNN_PARTICLES> & particles, const State * state ) {      
            if ( !state->recordViz.isEnviro ) {
                if ( getSampler() == NULL ) {
                    return 0;
                }
                return getSampler()->nnquery( hit, particles );
            }

            Vector3 dirFromEnviro = ( state->enviroDummySph.getCenter() - hit.position ).getNormalized();
            return getEnviroSampler()->nnquery( dirFromEnviro, particles );            
        }

        Particle getParticle( int index, const State * state ) const {            
            if ( state->recordViz.isEnviro ) {
                Particle p = getEnviroSampler()->getParticles()[ index ];                                
                //sampler->getSceneCenter() - enviroDummySph.getRadius() * recordPosition;
                Particle res = p;
                res.position = m_enviroSampler->getIncidentPosition( p );
                res.incidentDir = -p.position;
                res.distance *= m_enviroSampler->getDiscRadius();
                return res;
            }

            //IMPORTANCE_ASSERT( getSampler() != NULL );
            return getSampler()->particles[ index ];
        }

        /* It filters out particles that were really used for fitting */
        void getUsedParticles( const ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> & nnParticles, 
                int found, const Hit & hit, const State * state, IStack<Particle> & oUsedParticles ) const {
            if ( state->recordViz.isEnviro ) {
                for ( int i = 0; i < found; ++i ) {
                    const Particle p = getParticle( nnParticles[ i ].index, state );
                    oUsedParticles.push( p );
                }
                return;
            }

            /// We just need to fill in the oUsedParticles stack 
            class Wrapper {
            private:
                IStack<Particle> & m_data;
            public:
                Wrapper( IStack<Particle> & data ) : m_data( data ) {}
                void add( const Particle & p, const Vector3 & normal, Float ) {                              
                    Particle tmp = p;
                    // get back to incident direction
                    tmp.incidentDir = -tmp.incidentDir;
                    m_data.push( tmp );
                }
            } wrapper( oUsedParticles );
                        
            Config cfg;
            cfg.fitting.isFitOverSphere = false;
            DataFilter<Wrapper> filter( &m_sampler->particles, cfg );
            CacheStats dummy;
            filter.prepareData( nnParticles.ptr(), hit, found, dummy, wrapper );            
        }

        const Vector3 & getParticleRecVizDirection( int index, const State * state ) const {
            if ( state->recordViz.isEnviro ) {
                const Particle & p = getEnviroSampler()->getParticles()[ index ];
                return p.incidentDir;
            }

            return getSampler()->particles[ index ].incidentDir;
        }

    private:
        IEnviroSampler * m_enviroSampler;
        SimpleSampler<TDistributionModel> * m_sampler;
    };
}