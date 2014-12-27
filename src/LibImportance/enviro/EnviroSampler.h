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
#include "..\LibImportanceTypes.h"
#include "ImportanceSampler2d.h"
#include "..\caching\CachedSampler.h"
#include "ImportanceSampler1d.h"
#include "..\shared\buffers.h"
#include "ienviromap.h"
#include "pkddirections.h"
#include "parametricdirections.h" 
#include "..\caching\CacheStats.h"

namespace Importance {
    const bool  NON_UNIFORM_POSITIONS           = true;
    const float ONE_DEGREE_IN_EUCLIDIAN_DIST    = 0.017453f;

    /* Representation of importance function in the directional domain 
       (progressive kernel density or a parametric distribution) */
    typedef PKDDirections DirectionalImportance;
    //typedef ParametricDirections DirectionalImportance;

    //////////////////////////////////////////////////////////////////////////

    class IEnviroSampler {
    protected:
        BasicDistributionFactory * m_factory;
        // particle positions are DIRECTION TO ENVIRO
        IStack<Particle> m_particles;
    private:
        Vector3 m_sceneCenter;
        float m_discRadius;
    public:        
        virtual bool init( InputIterator& samples, const Config& config, Stats* stats, const EnviroMap* enviro ) = 0;

        virtual Vector3 chooseDir(const Vector2 random, float& pdf) const = 0;

        virtual float dirPdf(const Vector3 dirFromEnviro) const = 0;

        virtual Vector3 chooseOrigin(const Vector3 dirFromEnviro, const Vector2 random, float& pdf) = 0;

        virtual float originPdf(const Vector3 dirFromEnviro, const Vector3 origin) = 0;

        virtual void setDistributionFactory( BasicDistributionFactory * factory ) { m_factory = factory; }

        virtual std::string toString() const = 0;

        virtual VizAPI * getVizApi() = 0;

        virtual Vector3 getIncidentPosition( const Particle & p ) const = 0;

        virtual int nnquery( const Vector3 & dirFromEnviro, ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> & result ) const = 0;

        virtual bool isInitialized() const = 0;

        virtual void refreshSamples( InputIterator & samples, const EnviroMap * enviro = NULL, bool refineCache = true ) = 0;

        virtual const ImBitmap<float> * getDirectionsBmp() const = 0;
        virtual DebuggingBundle * getDebugMaps() = 0;

        virtual ~IEnviroSampler() {}

        const IStack<Particle> & getParticles() const {
            return m_particles;
        }

        IMPORTANCE_INLINE const Vector3 & getSceneCenter() const { return m_sceneCenter; }

        IMPORTANCE_INLINE float getDiscRadius() const { return m_discRadius; }

    protected:
        IMPORTANCE_INLINE void setSceneCenter( const Vector3 & center ) { m_sceneCenter = center; }

        IMPORTANCE_INLINE void setDiscRadius( float r ) { m_discRadius = r; }
    };

    //////////////////////////////////////////////////////////////////////////

    template <typename TDistributionModel>
    class EnviroSampler : public IEnviroSampler {
        template <class TDistributionModel> friend class CacheViz;
        template <class TDistributionModel> friend class RecordViz;
        template <class TDistributionModel> friend class SamplersFacade;
        friend struct State;
        template <class TDistributionModel> friend class EnviroSamplerVizAPI;
    protected:        

        // sampling directions FROM enviro
        class Sampler {
        public:
            ImportanceSampler1d<Float1Traits> columnSelector;
            IArray<ImportanceSampler1d<Float1Traits>> columns;
            float build(const ImBitmap<float>& bitmap, const Float1Traits& traits) {
                columns.resize(bitmap.getWidth());
                IArray<float> avgColumns(columns.size());
                for(int i = 0; i < columns.size(); ++i) {
                    auto it = bitmap.getDataIterator(i, 0);
                    const float avg = columns[i].construct(it, it+bitmap.getHeight(), traits);
                    avgColumns[i] = avg;
                }
                Float res = columnSelector.construct(avgColumns.begin(), avgColumns.end(), traits);
                IMPORTANCE_ASSERT( isReal( res ) );
                return res;
            }

            IMPORTANCE_INLINE Vector2 sample(const Vector2 random, float& pdf) const {
                Vector2 result;
                float dummy;
                result.x = columnSelector.sample(random.x, dummy);
                const int index = std::min(int(columns.size()-1), int(result.x*columns.size()));
                result.y = columns[index].sample(random.y, pdf);
                return result;
            }
        } directionSampler;

        IMPORTANCE_INLINE Hit getHit( const Vector3 & dirFromEnviro ) const {
            Hit res;
            res.position = dirFromEnviro;
            res.normal = Vector3( 0.f, 0.f, 1.f );
            return res;
        }

        virtual int nnquery( const Vector3 & dirFromEnviro, ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> & result ) const {
            Hit hit = getHit( dirFromEnviro );
            const float searchRadius = m_config.cache.maxEnviroRadius*ONE_DEGREE_IN_EUCLIDIAN_DIST;
            return m_tree.rangeQuery( hit.position, searchRadius, result.ptr(), m_config.particles.knn );
        }

        IMPORTANCE_INLINE Distribution* interpolate(const Vector3 dirFromEnviro, IResultBuffer& storage) {
            Hit dummyHit = getHit( dirFromEnviro );
            float dummy;
            Distribution* result = m_cache.getInterpolated(dummyHit, ONE_DEGREE_IN_EUCLIDIAN_DIST, m_config.cache.maxEnviroRadius, storage, dummy);
            
            if(result) {
                return result;
            } else {
                ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> particles;                
                const int found = nnquery( dirFromEnviro, particles );

                auto* rec = m_cache.createRecord(dummyHit);
                dummyHit.position = Vector3(NAN);
                dummyHit.normal = Vector3(0, 0, 1);
                CacheStats & cstats = rec->distr.m_cacheStats;
                cstats.particleSearchCount  = found;
                cstats.sqrSearchRadius      = particles[ 0 ].distSqr;
                if( TDistributionModel::getDistribution(m_factory, dummyHit, particles.ptr(), found, rec->distr) ) {
                    if ( rec->distr.m_cacheStats.particleSearchCount == 0 ) {
                        delete rec;
                        return NULL;
                    }
                    const TDistributionModel * dist = m_cache.add( rec, 
                                                                   ONE_DEGREE_IN_EUCLIDIAN_DIST,
                                                                   m_config.cache.maxEnviroRadius );
                    return new((void*)&storage) CrossfadeInterpolationResult<TDistributionModel>( dist );
                } else {
                    delete rec;
                    return NULL;
                }
            }
        }

        Vector3 _hemisphere2position(const Vector3 hemisphere, const Vector3 dirFromEnviro) const {
            const Frame rotator( -dirFromEnviro );
            /* If there is an invalid hemisphere position then generate ray position so that it possibly miss the scene */
            const Vector2 offset2 = hemisphere.z > 0.f ? hemisphereToDisk( hemisphere ) : Vector2( 10.f );
            //Vector3 offset = rotator.toWorld( Vector3(offset2.x, offset2.y, 0.f) * getDiscRadius() );
            //offset += getSceneCenter();
            Vector3 wOffset = rotator.toWorld( Vector3( offset2.x, offset2.y, 0.f ) );
            Vector3 offset = getSceneCenter() + ( wOffset - dirFromEnviro ) * getDiscRadius();  

//#ifdef CORONA_DEBUG            
#if 0
            if ( hemisphere.z > 0.f ) {
                const Vector3 reverse = _position2hemisphere( offset, dirFromEnviro );
                IMPORTANCE_ASSERT( (reverse.z == -1.f && hemisphere.z < 0.1e-6f) || 
                    abs( reverse.x-hemisphere.x ) < 5e-4f );
                IMPORTANCE_ASSERT( (reverse.z == -1.f && hemisphere.z < 0.1e-6f) ||
                    abs( reverse.y-hemisphere.y ) < 5e-4f );
                IMPORTANCE_ASSERT( (reverse.z == -1.f && hemisphere.z < 0.1e-6f) || 
                    abs( reverse.z-hemisphere.z ) < 5e-4f );
            }
#endif
            IMPORTANCE_ASSERT( offset.isReal() );
            return offset;
        }

        Vector3 _position2hemisphere(const Vector3 position, const Vector3 dirFromEnviro) const {
            const Frame rotator(-dirFromEnviro);
            Vector3 offset = (position - getSceneCenter());
            offset = rotator.toLocal(offset);
            Vector2 diskOffset(offset.x, offset.y);
            diskOffset /= getDiscRadius();
            if ( diskOffset.square() > 1.f ) {
                // position is outside the disk
                return Vector3( 0.f, 0.f, -1.f );
            }
            const Vector3 result = diskToHemisphere(diskOffset);
            return result;
        }

        void fill( IStack<Particle> & particles, InputIterator & samples ) {
            while( samples.isValid() ) {
                Particle p;
                p.weight    = samples.weight();        
                p.normal    = Vector3(0, 0, 1);           
                p.incidentDir = -this->_position2hemisphere(samples.position(), -samples.incidentDirection());
                // if this assert fails it probably means, that disk radius does not cover the whole scene
                IMPORTANCE_ASSERT( p.incidentDir.x != 0.f || p.incidentDir.y != 0.f || p.incidentDir.z != 1.f );
                p.position = -samples.incidentDirection();

                /* We compute distance of particle here. We ignore passed samples.distance() */
                Vector3 fromPosition = samples.position();
                Vector3 toPosition   = getIncidentPosition( p );
                p.distance           = ( toPosition - fromPosition ).length() / getDiscRadius();
                IMPORTANCE_ASSERT( isReal( p.distance ) );

                particles.push( p );
                samples.next();                
            }
        }

        IMPORTANCE_INLINE void resetDirectionSampler( ImBitmap<float> rasterPdf ) {
            this->directionSampler = EnviroSampler::Sampler();
            this->m_avgDirValue    = this->directionSampler.build( rasterPdf, Float1Traits());            
        }

    public:

        EnviroSampler() : initialized(false), m_directionalImportance( m_particles, m_tree ) {}

        virtual ~EnviroSampler() {}

        virtual bool isInitialized() const {
            return this->initialized;
        }

        /* Reconstructs original position of a particle (up to some error because now the position if not on the environment sphere but on the tangent disc )
            We need this because of raping original data items in Particle structure. Particle::position is
            used for storing direction "fromEnviroMap", i.e. opposite direction to particle's incident direction
         */
        virtual Vector3 getIncidentPosition( const Particle & p ) const {
            return this->_hemisphere2position( -p.incidentDir, p.position );
        }

        // inputIterator jsou fotony, tedy smeruji DO enviro
        virtual bool init( InputIterator& samples, const Config& config, Stats* stats, const EnviroMap* enviro ) {
            setSceneCenter( enviro->getSceneCenter() );
            setDiscRadius( enviro->getDiscRadius() );            
            this->m_stats = stats;
            this->m_config = config;
            this->m_config.misCorrection.use = false;
            this->m_config.perfectSampling = false;
            this->m_config.sceneBboxMin = Vector3(-1.f);
            this->m_config.sceneBboxMax = Vector3(1.f);            
                
            fill( m_particles, samples );

            if(m_particles.size() == 0) {
                return false; 
            }
            this->m_tree.build(this->m_particles.begin(), this->m_particles.end(), PARTICLE_TREE_LEAF_SIZE);
    
            m_enviroMap.resize( 1000, 500 );
#ifdef LIBIMP_STATS
            m_debugMaps.resize( 1000, 500 );
#endif
            m_directionalImportance.init( enviro, m_enviroMap, m_debugMaps, m_config );
            resetDirectionSampler( m_enviroMap );

            this->m_factory->init(this->m_particles, &this->m_config);
            //const BoundingBox3 bbox = BoundingBox3(config.sceneBboxMin, config.sceneBboxMax);
            BoundingBox3 bbox;
            for( auto it = m_particles.begin(); it < m_particles.end(); ++it ) { bbox += it->position; }
            this->m_cache.init(this->m_config, stats, bbox);
            this->initialized = true;            

            return true;
        }

        virtual Vector3 chooseDir(const Vector2 random, float& pdf) const {
            IMPORTANCE_ASSERT(isInitialized());

            const Vector2 tmp = this->directionSampler.sample(random, pdf);
            const Vector3 result = EnviroSamplerMapping::getDir(tmp);
            pdf = pdf/(this->m_avgDirValue*4*PI);

            //const Vector3 result = uniformSphere( random.x, random.y );
            //pdf = dirPdf( result );
            IMPORTANCE_ASSERT(result.isNormalized());
            return result;
        }

        virtual float dirPdf(const Vector3 dirFromEnviro) const {
            IMPORTANCE_ASSERT(isInitialized());

            //return 1.f / ( 4.f * PI );

            const Vector2 randCoords = EnviroSamplerMapping::getDirInverse(dirFromEnviro);
            const int x = std::min(this->m_enviroMap.getWidth()-1, int(this->m_enviroMap.getWidth()*randCoords.x));
            const int y = std::min(this->m_enviroMap.getHeight()-1, int(this->m_enviroMap.getHeight()*randCoords.y));
            const float value = this->m_enviroMap(x, y);
            IMPORTANCE_ASSERT(isReal(value/m_avgDirValue) && value/m_avgDirValue >= 0.f);
            return value/(4*PI*this->m_avgDirValue);
        }

        virtual Vector3 chooseOrigin(const Vector3 dirFromEnviro, const Vector2 random, float& pdf) {
            IMPORTANCE_ASSERT(isInitialized());
            ResultBuffer resultBuffer;
            Distribution* distr = this->interpolate(dirFromEnviro, resultBuffer);
            Vector3 result;
            if(distr && NON_UNIFORM_POSITIONS) {
                float cachePdf;
                Vector3 dir;
                distr->sampleDirections( &random, &dir, &cachePdf, 1 );
                result = _hemisphere2position( dir, dirFromEnviro );
                pdf = cachePdf * (2*PI) / (PI * getDiscRadius() * getDiscRadius());
                IMPORTANCE_ASSERT( result.isReal() );                
                IMPORTANCE_ASSERT( pdf > 0.f || ( pdf == 0.f && ( (getSceneCenter() - dirFromEnviro - result ).size() / getDiscRadius()) > 1.f ) );
            } else {
                result = _hemisphere2position( uniformHemisphereShirley(random.x, random.y), dirFromEnviro );             
                pdf = 1/(PI * getDiscRadius() * getDiscRadius());
            }
            IMPORTANCE_ASSERT( pdf >= 0.f );
            return result;
        }

        virtual float originPdf(const Vector3 dirFromEnviro, const Vector3 origin) {
            IMPORTANCE_ASSERT(isInitialized());
            ResultBuffer resultBuffer;
            Distribution* distr = this->interpolate(dirFromEnviro, resultBuffer);
            if(distr && NON_UNIFORM_POSITIONS) {
                const Vector3 onHemisphere = this->_position2hemisphere(origin, dirFromEnviro);
                IMPORTANCE_ASSERT(onHemisphere.isNormalized());
                float cachePdf;
                distr->pdfs(&onHemisphere, &cachePdf, 1);
                IMPORTANCE_ASSERT( isReal( cachePdf ) && !isNaN( cachePdf ) );
                IMPORTANCE_ASSERT( onHemisphere.z >= 0.f || cachePdf == 0.f );
                return cachePdf*2*PI/(PI * getDiscRadius() * getDiscRadius());
            } else {
                return 1/(PI * getDiscRadius() * getDiscRadius());        
            }
        }    

        virtual std::string toString() const {
            std::ostringstream ostr;
            ostr << "EnviroSamler [" << std::endl;
            ostr << "Number of particles: " << m_particles.size() << std::endl;
            ostr << "Is initialized: " << ( initialized ? "true" : "false" ) << std::endl;
            ostr << "Disk radius: " << getDiscRadius() << std::endl;
            ostr << "Scene center: " << getSceneCenter().toString() << std::endl;
            ostr << "Average direction value: " << m_avgDirValue << std::endl;
            ostr << "------------------------------" << std::endl;
            ostr << m_directionalImportance.toString() << std::endl;
            ostr << "]" << std::endl;
            return ostr.str();          
        }



        virtual void refreshSamples( InputIterator & samples, const EnviroMap * enviro = NULL, bool refineCache = true ) {            
            IMPORTANCE_ASSERT( this->m_config.particles.knn <= MAX_KNN_PARTICLES );

            if ( !isInitialized() || !samples.isValid() ) {
                return;
            }
             
            if ( m_config.fitting.enviroType == EStatic ) {
                ILog( EWarn, "Progressive update of both directions and positions on the environment map is switched off!" );
                return;
            }            
            
            /* Store new particles and build tree */
            m_particles.clear();
            fill( m_particles, samples );            
            m_tree = PointKdTree<Kd3PositionTraits>();
            m_tree.build( m_particles.begin(), m_particles.end(), PARTICLE_TREE_LEAF_SIZE );

            /* Importance function in a directional domain, multiply with the environment map */
            if ( m_config.fitting.enviroType == EDynamic || m_config.fitting.enviroType == EStaticPositions ) {
                m_directionalImportance.refresh( enviro, m_enviroMap, m_debugMaps );
                resetDirectionSampler( m_enviroMap );
                ILog( EInfo, "Importance envmap: directions density updated." );
            }            

            /* Importance function of positions */
            if ( m_config.fitting.enviroType == EDynamic || m_config.fitting.enviroType == EStaticDirections ) {
                int counter = 0;
                auto & recs = m_cache.getRecords();

                #pragma omp parallel for
                for ( int i = 0; i < (int) recs.size(); ++i ) {
                    CacheRecord<TDistributionModel> * rec = recs[ i ];
                    counter++;
                    ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> searchResults;
                    Hit hit;
                    hit.normal      = rec->normal;
                    hit.position    = rec->position;
                    const int nFound = nnquery( hit.position, searchResults );                    
                    CacheStats & cstats = rec->distr.m_cacheStats;
                    cstats.particleSearchCount  = nFound;
                    cstats.sqrSearchRadius      = searchResults[ 0 ].distSqr;
                    TDistributionModel::getDitributionUpdated( m_factory, hit, searchResults.ptr(), nFound, rec->distr);
                    if ( refineCache /*&& nFound > 0*/ ) {    
                        m_cache.updateRecordRadius( rec, ONE_DEGREE_IN_EUCLIDIAN_DIST, m_config.cache.maxEnviroRadius );
                    }
                }
                if ( refineCache ) {
                    m_cache.rebuildTree();
                }
                ILog( EInfo, "Importance envmap: cache of position records updated - approx. %d records were updated", (int) counter );
            }            
        }

        virtual VizAPI * getVizApi() { return new EnviroSamplerVizAPI<TDistributionModel>( this ); }

        virtual const ImBitmap<float> * getDirectionsBmp() const {            
            return &m_enviroMap; 
        }
        virtual DebuggingBundle * getDebugMaps() {             
            return &m_debugMaps; 
        }
    
        //////////////////////////////////////////////////////////////////////////
        // Data
        //////////////////////////////////////////////////////////////////////////

private:
        Stats* m_stats;
        Config m_config;

        ImBitmap<float> m_enviroMap;
        DebuggingBundle m_debugMaps;
        float m_avgDirValue;        

        // vsechno v cache je ulozeny s position = dirfromenviro
        ImportanceCache<TDistributionModel> m_cache;
        PointKdTree<Kd3PositionTraits> m_tree;
        DirectionalImportance m_directionalImportance;
        bool initialized;        
    };
}