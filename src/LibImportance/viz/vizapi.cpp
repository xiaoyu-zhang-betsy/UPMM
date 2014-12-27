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


#define __WXMSW__
#define WXUSINGDLL
#include <wx/wx.h>
#include <gl/GLU.h>

#include "vizapi.h"
#include "../caching/CachedSampler.h"
#include "../enviro/EnviroSampler.h"
#include "../viz/State.h"

#include "../distributionmodels.h"


namespace Importance {

    template<class TDistributionModel>
    VizHitRecord EnviroSamplerVizAPI<TDistributionModel>::getVizHitRecord( const Hit & hit, IResultBuffer & buffer, const State & state )
    {
        VizHitRecord result;
        if ( m_sampler == NULL || !m_sampler->isInitialized() ) {
            result.isValid = false;
            return result;
        }

        result.selectedId = -1;        
        float d = INFINITY;
        m_sampler->m_cache.lock.lockRead();        
        for(int i = 0; i < m_sampler->m_cache.records.size(); ++i) {
            const float newD = (m_sampler->m_cache.records[i]->position-hit.normal).size();
            if(newD < d) {
                result.selectedId = i;
                d = newD;
            }
        }
        m_sampler->m_cache.lock.unlockRead();

        if ( result.selectedId == -1 ) {
            result.isValid = false;
            return result;
        }

        const auto * cacheRecord = m_sampler->m_cache.records[ result.selectedId ];
        result.isValid           = true;
        result.pDistr            = new (&buffer) ImmediateDistribution<TDistributionModel>( cacheRecord->distr );
        result.normal            = cacheRecord->normal;
        result.position          = recordToPosition( state, cacheRecord );
        return result;
    }


    template<class TDistributionModel>
    VizHitRecord CachedSamplerVizAPI<TDistributionModel>::getVizHitRecord( const Hit & hit, IResultBuffer & buffer, const State & state )
    {
        VizHitRecord result;           
        if ( m_sampler == NULL ) {
            result.isValid = false;
            return result;
        }

        result.selectedId = -1;
        
        float d = INFINITY;
        if ( !state.interpolateOnClick ) {
            m_sampler->cache.lock.lockRead();        
            for(int i = 0; i < m_sampler->cache.records.size(); ++i) {
                const auto rec = m_sampler->cache.records[ i ];
                const float newD = (rec->position - hit.position).size();
                if(newD < d && dot( rec->normal, hit.normal ) > 0.f ) {
                    result.selectedId = i;
                    d = newD;
                }
            }
            m_sampler->cache.lock.unlockRead();

            if ( result.selectedId != -1 ) {
                const auto * cacheRecord  = m_sampler->cache.records[ result.selectedId ];
                result.pDistr             = new (&buffer) ImmediateDistribution<TDistributionModel>( cacheRecord->distr ) ;
                result.normal             = cacheRecord->normal;
                result.position           = cacheRecord->position;
                if ( result.selectedId == -1 ) {                
                    result.isValid = false;
                    return result;
                }
                result.isValid = true;
            }
        }
        
        if ( result.selectedId == -1 ) {
            /* If interpolateOnClick or there are no records then try to get it from sampler 
               near the hit point */
            result.pDistr   = m_sampler->getDistribution( hit, buffer );
            result.isValid  = result.pDistr != NULL;
            result.normal   = hit.normal;
            result.position = hit.position;
        }        

        return result;
    }

    template<class TDistributionModel>
    VizHitRecord SimpleSamplerVizAPI<TDistributionModel>::getVizHitRecord( const Hit & hit, IResultBuffer & buffer, const State & state )
    {
        VizHitRecord result;
        if ( m_sampler == NULL ) {
            result.isValid = false;
            return result;
        }

        result.isValid     = true;
        result.pDistr      = m_sampler->getDistribution( hit, buffer );
        result.normal      = hit.normal;
        result.position    = hit.position;
        return result;
    }

    template<class TDistributionModel>
    SimpleSamplerVizAPI<TDistributionModel>::SimpleSamplerVizAPI( SimpleSampler<TDistributionModel> * sampler ) : m_sampler( sampler )
    {

    }

    template<class TDistributionModel>
    CachedSamplerVizAPI<TDistributionModel>::CachedSamplerVizAPI( CachedSampler<TDistributionModel> * sampler ) : m_sampler( sampler )
    {

    }

    template<class TDistributionModel>
    EnviroSamplerVizAPI<TDistributionModel>::EnviroSamplerVizAPI( EnviroSampler<TDistributionModel> * sampler ) : m_sampler( sampler )
    {

    }


    template<class TDistributionModel>
    void CachedSamplerVizAPI<TDistributionModel>::draw( const State & state ) const
    {
        if ( state.showRecords ) {
            paintAllNonEnviroRecords( state );
        }

        if ( state.interpolateOnClick && state.recordViz.vRec.isValid ) {
            paintActiveRecords( state );
        }
    }

    template<class TDistributionModel>
    void CachedSamplerVizAPI<TDistributionModel>::paintActiveRecords( const State & state ) const 
    {
        if ( m_sampler == NULL || !state.recordViz.vRec.isValid || state.recordViz.isEnviro ) {
            return;
        }

        glColor3f(0.f, 1.f, 0.f);
        glLineWidth(2.f);
        glBegin(GL_LINES);
         InterpolationResult<TDistributionModel> * clickedDistribution = 
            dynamic_cast<InterpolationResult<TDistributionModel>*>( state.recordViz.vRec.pDistr );
        
        for(int i = 0; i < m_sampler->cache.records.size(); ++i) {
            const CacheRecord<TDistributionModel>& actual = *m_sampler->cache.records[i];
            if( clickedDistribution ) {
                const int recordsUsed = clickedDistribution->numRecordsUsed();
                for(int j = 0; j < recordsUsed; ++j) {
                    if(&actual.distr == clickedDistribution->getModel(j)) {
                        glVertex3f(state.lastClick.position.x, state.lastClick.position.y, state.lastClick.position.z);
                        glVertex3f(actual.position.x, actual.position.y, actual.position.z);
                        glEnd();
                        VizAPI::drawRadius(actual.position, actual.normal, actual.validRadius());
                        VizAPI::drawRadius(actual.position, actual.normal, actual.originalRadius, Vector3(0,1,0));
                        glBegin(GL_LINES);
                        glColor3f(0.f, 1.f, 0.f);
                    }
                }
            }
        }
        glEnd();
    }


    void VizAPI::drawRadius( const Vector3 position, const Vector3 normal, const float radius, const Vector3 c /*= Vector3(1,0,0) */ )
    {
        Frame rotator(normal);

        glColor3f(c.x, c.y, c.z);

        const int SEGMENTS = 32;
        glBegin(GL_LINE_LOOP);
        for(int i = 0; i < SEGMENTS; ++i) {
            const float angle = 2*PI*i/float(SEGMENTS);
            Vector3 pos = rotator.toWorld(Vector3(cos(angle), sin(angle), 0.001f))*radius + position;
            glVertex3f(pos.x, pos.y, pos.z);
        }
        glEnd();
    }


    template<class TDistributionModel>
    void CachedSamplerVizAPI<TDistributionModel>::paintAllNonEnviroRecords( const State & state ) const
    {
        if ( m_sampler == NULL ) {
            return;
        }

        glPointSize(2.0f);
        glBegin(GL_POINTS);
        for(int i = 0; i < m_sampler->cache.records.size(); i += state.everyNthRecord) {
            const auto& record = m_sampler->cache.records[ i ];
            
            setRecordColor(*record, state);
            glVertex3f(record->position.x, record->position.y, record->position.z);

            if( !state.interpolateOnClick && i == state.recordViz.vRec.selectedId 
                && !state.recordViz.isEnviro ) {
                glEnd();
                VizAPI::drawRadius(record->position, record->normal, record->validRadius());                    
                VizAPI::drawRadius(record->position, record->normal, record->originalRadius, Vector3(0,1,0));         
                glBegin(GL_POINTS);
            }
        }
        glEnd();

        if(state.showRadii) {
            for(int i = 0; i < m_sampler->cache.records.size(); i += state.everyNthRecord) {
                const auto& record = m_sampler->cache.records[ i ];            
                VizAPI::drawRadius(record->position, record->normal, record->validRadius());           
            }
        }
        
        if(state.isDisplayNormals) {
            glBegin(GL_LINES);
            for(int i = 0; i < m_sampler->cache.records.size(); i += state.everyNthRecord) {
                const auto& record = m_sampler->cache.records[ i ];
                setRecordColor(*record, state);
                Vector3 to = record->position + 0.01f * state.bbox.size() * record->normal;
                Vector3 from = record->position;
                glVertex3f(from.x, from.y, from.z);
                glVertex3f(to.x, to.y, to.z);
            }
            glEnd();
        }
    }


    template<class TDistributionModel>
    void EnviroSamplerVizAPI<TDistributionModel>::paintAllEnviroRecords( const State & state ) const
    {
        if ( m_sampler == NULL || !m_sampler->isInitialized() ) {
            return;
        }

        glPointSize( 2.0f );
        glBegin( GL_POINTS );
        for ( int i = 0; i < m_sampler->m_cache.records.size(); i += state.everyNthRecord ) {
            const auto& record = m_sampler->m_cache.records[i];
            Vector3 normal = state.lastClick.normal;
            Vector3 pos = recordToPosition( state, record );

            setRecordColor(*record, state);

            if(!state.interpolateOnClick && i == state.recordViz.vRec.selectedId
                && state.recordViz.isEnviro ) {
                glEnd();                
                VizAPI::drawRadius(pos, normal, m_sampler->getDiscRadius() * record->validRadius() );
                VizAPI::drawRadius(pos, normal, m_sampler->getDiscRadius() * record->originalRadius, Vector3(0,1,0));
                VizAPI::drawRadius(pos, normal, m_sampler->getDiscRadius(), Vector3(0,1,0));
                glBegin(GL_POINTS);
            }
            glVertex3f(pos.x, pos.y, pos.z);
        }
        glEnd();

        if ( state.showRadii ) {
            for ( int i = 0; i < m_sampler->m_cache.records.size(); i += state.everyNthRecord ) {
                const auto& record = m_sampler->m_cache.records[i];                                    
                Vector3 normal = state.lastClick.normal;
                Vector3 pos = state.enviroDummySph.getCenter() - state.enviroDummySph.getRadius() * record->position;
                VizAPI::drawRadius(pos, normal, record->validRadius());                                
            }
        }
    }


    template<class TDistributionModel>
    void EnviroSamplerVizAPI<TDistributionModel>::draw( const State & state ) const
    {
        if ( state.showRecords ) {
            paintAllEnviroRecords( state );
        }
    }

    template <class TDistributionModel>
    Vector3 EnviroSamplerVizAPI<TDistributionModel>::recordToPosition( const State & state, 
        const CacheRecord<TDistributionModel> * rec ) {
        return state.enviroDummySph.getCenter() - state.enviroDummySph.getRadius() * rec->position;
    }

    template class SimpleSamplerVizAPI<Jensen>;
    template class SimpleSamplerVizAPI<Hey>;
    template class SimpleSamplerVizAPI<Pharr>;
    template class SimpleSamplerVizAPI<GaussianMixtureSSEHemisphere>;

    template class CachedSamplerVizAPI<Jensen>;
    template class CachedSamplerVizAPI<Hey>;
    template class CachedSamplerVizAPI<Pharr>;
    template class CachedSamplerVizAPI<GaussianMixtureSSEHemisphere>;

    template class EnviroSamplerVizAPI<Jensen>;
    template class EnviroSamplerVizAPI<Hey>;
    template class EnviroSamplerVizAPI<Pharr>;
    template class EnviroSamplerVizAPI<GaussianMixtureSSEHemisphere>;
}