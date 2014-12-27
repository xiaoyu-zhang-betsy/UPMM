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

#ifndef ___CACHED_SAMPLER___
#define ___CACHED_SAMPLER___

#include "DynamicKdTree.h"
#include "../shared/StaticArray.h"
#include "../shared/Vector2.h"
#include "../shared/Utils.h"
#include "ReaderWriterLock.h"
#include "../simple/SimpleSampler.h"
#include "MultirefOctree.h"
#include "../shared/IIterator.h"
#include "CacheTypes.h"

#pragma warning(push)
#pragma warning(disable:4127)


namespace Importance {

#define USE_OCTREE

    enum InterpolationMetricType {
        METRIC_KL_DIVERGENCE = 0,
        METRIC_IRRADIANCE = 1,
        METRIC_D_IRRADIANCE = 2,
        METRIC_NONE = 666,
        METRIC_TEST = 5,
        METRIC_GRADIENT = 6
    };

    const int NEIGHBOUR_CLAMPING_KNN = 20;

    template<class TDistributionModel>
    class ImportanceCache {
        friend void testViz();
        template <class TDistributionModel> friend class CacheViz;
        template <class TDistributionModel> friend class RecordViz;
        friend struct State;        
    public:
        typedef CacheRecord<TDistributionModel> Record;       

        mutable ReaderWriterLock lock;

        // all records in simple linear array
        IStack<Record*> records;
    protected:

        const Config* config;

#ifdef USE_OCTREE
        MultirefOctree octree;
        typedef MultirefOctree::IIterator FoundRecordsIterator;
#else
        ImDynamicKdTree<DynamicKd3PtrPositionTraits> tree;
        int newlyInserted;
        static const int LEAF_SIZE = 25;
        class FoundRecordsIterator {
        public:
            FoundRecordsIterator( const ImStaticArray<KdQueryResult, NEIGHBOUR_CLAMPING_KNN> & nearest, int found ) : m_nearest( nearest ), m_found( found ), m_position( 0 ) {}
            IMPORTANCE_INLINE int getNext() { 
                IMPORTANCE_ASSERT( m_position < m_found );
                return m_nearest[ m_position ].index; 
            }
            IMPORTANCE_INLINE bool hasNext() const { return m_position < m_found; }
        private:
            void operator=( const FoundRecordsIterator& ) {}
            const ImStaticArray<KdQueryResult, NEIGHBOUR_CLAMPING_KNN> & m_nearest;
            int m_position;
            int m_found;
        };
#endif

        Stats* stats;
    protected:   
        struct RecordSearchResult {                
            int index;
            float distances;
            float error;
            bool rejected;
            IMPORTANCE_INLINE bool operator<(const RecordSearchResult& other) const {
                return error > other.error;
            }
        }; 

        struct ValidityCriteria {        
            Float R_alpha;
            Float R_min_photon;
            Float R_max_photon;
            Float R_near;
            Float R_tentative;
            Float minValidity;
            Float maxValidity;

#ifdef LIBIMP_DEBUG
            ValidityCriteria() {
                R_alpha = R_min_photon = R_max_photon = R_near = R_tentative = minValidity = maxValidity = NAN;
            }
#endif

            IMPORTANCE_INLINE Float evaluateValidityRadius( const Record * record ) const {
                IMPORTANCE_ASSERT( this->isValid() );
                const Float newRadius = std::min( std::max(R_min_photon, R_tentative), record->validRadius() );
                const Float res = clamp(newRadius, minValidity, maxValidity);
                IMPORTANCE_ASSERT( res >= 0.f );
                return res;
            }

            IMPORTANCE_INLINE bool isValid() const {
                return !isNaN( R_alpha ) &&
                       !isINF( R_alpha ) &&
                       !isNaN( R_min_photon ) &&
                       !isINF( R_min_photon ) &&
                       !isNaN( R_max_photon ) &&
                       !isINF( R_max_photon ) &&                                              
                       
                       !isNaN( R_near ) &&

                       !isNaN( R_tentative ) &&
                       !isINF( R_tentative ) &&
                       !isNaN( minValidity ) &&
                       !isINF( minValidity ) &&
                       !isNaN( maxValidity ) &&
                       !isINF( maxValidity );
            }
        };

        IMPORTANCE_INLINE void computeValidityCriteria( Float maxRadius, Float pixel2world, 
            bool changeOriginalRadius, Record * record, ValidityCriteria & c ) {
            const CacheStats & cstats = record->distr.m_cacheStats;
            IMPORTANCE_ASSERT( isReal( cstats.sqrSearchRadius ) && cstats.particleSearchCount >= 0 );
            Float photonSearchRadius = std::sqrt( cstats.sqrSearchRadius );

            c.R_alpha       = TDistributionModel::validityRadius(record->distr, *config);
            c.R_min_photon  = config->cache.photonClampingMin*photonSearchRadius;
            c.R_max_photon  = config->cache.photonClampingMax*photonSearchRadius;
            c.minValidity   = pixel2world*this->config->cache.minRadius;
            c.maxValidity   = pixel2world*maxRadius;
            /* Generalization achieved by parametric representation cannot be trusted yet because we have observed 
               very low number of samples. Thus restrict the radius by particle search radius. */
            if ( record->distr.isUndersampled() ) {
                c.R_min_photon = photonSearchRadius;
            }
            /* There is no reliable information that could be used to set the validity radius (not even to use particle 
               search radius). Thus we use an user 
               specified threshold based on the scene size. */
            if ( cstats.particleSearchCount < 5 ) {
                c.R_min_photon = c.maxValidity;
            }

            IMPORTANCE_ASSERT(isReal(c.R_alpha) && record->normal.isReal());
            
            c.R_near = INFINITY;

            if( config->cache.neighbourClamping ) {                                
                lock.lockRead();
#ifdef USE_OCTREE
                FoundRecordsIterator it( &this->octree, record->position );
                FoundRecordsIterator it2 = it;
#else
                ImStaticArray<KdQueryResult, NEIGHBOUR_CLAMPING_KNN> nearest;
                const int found = this->tree.getClosest(record->position, INFINITY, NEIGHBOUR_CLAMPING_KNN, nearest.ptr());
                FoundRecordsIterator it( nearest, found );
                FoundRecordsIterator it2( nearest, found );
#endif
                c.R_near = clampRecordByNeighbours( it, record );
                c.R_tentative = std::min(c.R_near, std::min(c.R_max_photon, c.R_alpha));
                // this is split in two phases so that we already know all 
                clampNeighboursByRecord( it2, c, record, changeOriginalRadius );
                lock.unlockRead();
            }

            c.R_tentative = std::min(c.R_near, std::min(c.R_max_photon, c.R_alpha));
        }

        IMPORTANCE_INLINE void takeTheBestResult( const int recordIndex, const Float maxDistanceInverse, const Hit & queryHit, RecordSearchResult & res ) const {
            const Record * record       = records[ recordIndex ];
            float distance              = (record->position - queryHit.position).sizeApprox();
            const float distanceTerm    = distance * maxDistanceInverse;
            const float normalTerm      = 2 * Ff::sqrtFast(1.001f - dot(record->normal, queryHit.normal));
            const float error           = distanceTerm + normalTerm;
            IMPORTANCE_ASSERT(isReal(error) && error >= 0);

            if ( error < res.error ) {
                res.index       = recordIndex;
                res.distances   = distance;
                res.error       = error;
                res.rejected    = distance > record->validRadius() || normalTerm > 1.f;                
            }
        }

    public:

        const IStack<Record*>& getRecords() const {
            return records;
        }

        /// developement hack
        const IStack<Record*>& getRecords() {
            return records;
        }

        virtual ~ImportanceCache() {
            for(auto it = records.begin(); it != records.end(); it++) {
                delete *it;
            }
        }

        // initializes the cache for given scene
        void init(const Config& config, Stats* stats, const BoundingBox3 & bbox) {
            this->stats = stats;
            this->config = &config;            
#ifdef USE_OCTREE
            this->octree.build(this->records, bbox, 1/this->config->cache.invGeometryError);
#else 
            this->tree.build(this->records.cbegin(), this->records.cend(), bbox, LEAF_SIZE);
#endif
        }

        void fillRecords(ImStaticArray<const TDistributionModel*, MAX_CACHE_KNN> & distModels, const KdQueryResult * results, const int used) const {
            for(int i = 0; i < used; ++i) {
                distModels[i] = &records[results[i].index]->distr;
            }
        }


        IMPORTANCE_INLINE Record * createRecord( const Hit & hit ) const {
            Record * res = new Record(Vector3(hit.position), Vector3(hit.normal), -1.f);
            res->setValidRadius( INFINITY );
            return res;
        }


        const TDistributionModel* add(Record * record, const Float pixel2world, 
const Float maxRadius ) {
                IMPORTANCE_ASSERT( record->validRadius() == INFINITY );
                ValidityCriteria c;
                computeValidityCriteria( maxRadius, pixel2world, true, record, c );                
                record->setValidRadius( c.evaluateValidityRadius( record ) );
				record->originalRadius = record->validRadius();

                this->lock.lockWrite();
                IMPORTANCE_INC_STATS(stats->cache.records);
                IMPORTANCE_INC_STATS(stats->cache.newRecords);
                if ( record->distr.isUndersampled() ) {
                    IMPORTANCE_INC_STATS(stats->cache.newUndersampledRecords);
                }
                this->records.push(record);
#ifdef USE_OCTREE
                this->octree.addRecord(record->position, record->validRadius(), 1/this->config->cache.invGeometryError);
#else
                this->tree.add(record->position, int(this->records.size()-1));
                ++this->newlyInserted;
                if(newlyInserted*2 > this->records.size() && this->records.size() > 1000) {
                    BoundingBox3 bbox(this->config->sceneBboxMin, this->config->sceneBboxMax);
                    this->tree.build(this->records.cbegin(), this->records.cend(), bbox, LEAF_SIZE);
                    this->newlyInserted = 0;
#ifdef LIBIMP_DEBUG
                    std::cout << "---------Rebuilding tree\n";
#endif
                }
#endif
                this->lock.unlockWrite();
                return &record->distr;
        }



		// This method serves for record validity radius reduction 
		// Does not use finiteGradient -> hard to implement and not used anymore
		void updateRecordRadius( Record * record, Float pixel2world, const Float maxRadius ) {
            record->distr.m_cacheStats.sqrSearchRadius = record->distr.getPhotonMaxDistance();
            ValidityCriteria c;
            computeValidityCriteria( maxRadius, pixel2world, false, record, c );            
            record->setValidRadius( c.evaluateValidityRadius( record ) );			
		}
      

        InterpolationResult<TDistributionModel>* getInterpolated(const Hit& hit, const float pixel2world,
            const Float maxRadius, IResultBuffer& buffer, float& outMetric) const {

            const float maximalDistance = pixel2world * maxRadius;
            const float maxDistanceInverse = 1/maximalDistance;

            RecordSearchResult res;
            res.error = INFINITY;
                    
            this->lock.lockRead();
#ifdef USE_OCTREE
            MultirefOctree::IIterator it(&this->octree, hit.position);
            int count = 0;
            while(it.hasNext()) {
                count++;                
                takeTheBestResult( it.getNext(), maxDistanceInverse, hit, res );
            }
            const int found = count;
            IMPORTANCE_INC_STATS_MORE(stats->cache.recordsFound, count);
#else 
            const int KNN_SEARCHED = 8;
            ImStaticArray<KdQueryResult, KNN_SEARCHED> knnResults;
            const int found = tree.getClosest(hit.position, maximalDistance, KNN_SEARCHED, knnResults.ptr());
            for(int i = 0; i < found; ++i) {                
                takeTheBestResult( knnResults[i].index, maxDistanceInverse, hit, res );
            }
            IMPORTANCE_INC_STATS_MORE(stats->cache.recordsFound, found);
#endif

            if(found == 0 || res.rejected) {                
                lock.unlockRead();
                return NULL;
            }

            const Record* recordA = records[res.index];
            lock.unlockRead();

//            if( res.distances*2 < recordA->validRadius() ) {
//                return new((void*)&buffer) CrossfadeInterpolationResult<TDistributionModel>(&recordA->distr) ;
//            }
            
//            return NULL;

              return new((void*)&buffer) CrossfadeInterpolationResult<TDistributionModel>(&recordA->distr) ;
        }


        void rebuildTree() {
            const BoundingBox3 bbox = BoundingBox3(this->config->sceneBboxMin, this->config->sceneBboxMax);
#ifdef USE_OCTREE
            this->octree.build(this->records, bbox, 1/this->config->cache.invGeometryError);
#else 
            newlyInserted = 0;
            this->tree.build(this->records.cbegin(), this->records.cend(), bbox, LEAF_SIZE);
#endif
        }

    protected:

        Float clampRecordByNeighbours( FoundRecordsIterator & nearest, const Record * newRecord ) const {
            Float R_near = INFINITY;            

            /* Check neighbours and clamp radius of newly added record so that at maximum it reaches to the closest further edge of a record. */                                    
            while( nearest.hasNext() ) {
                const Record* actual = this->records[ nearest.getNext() ];
                if(dot(actual->normal, newRecord->normal) > 0.4f) {
                    R_near = std::min(R_near, (actual->position - newRecord->position).sizeApprox() + actual->validRadius());
                }
            }            
            return R_near;
        }

        /* changeOriginalRadius has meaning only for debugging in the visualization tool. If nearest records are clamped 
           and changeOriginalRadius is true than we set this clamped radius as the original radius. Thus this should be 
           false if neighbourClamping is used in radius update. */
        void clampNeighboursByRecord( FoundRecordsIterator & nearest, const ValidityCriteria & c, 
            const Record * record, bool changeOriginalRadius ) {
            IMPORTANCE_ASSERT( c.isValid() );            

            const Float actualMin = c.evaluateValidityRadius( record );

            /* Now check all nearest records if they need to be clamped due to adding the new record. */
            while( nearest.hasNext() ) {
                Record* actual = this->records[ nearest.getNext() ];
                const float upperBound = actualMin + (actual->position - record->position).sizeApprox();
                if(dot(actual->normal, record->normal) > 0.4f && actual->validRadius() > upperBound) {
                    actual->setValidRadius( upperBound );
                    if ( changeOriginalRadius ) {
                        actual->originalRadius = upperBound;
                    }
                }
            }                                
        }
        
    };


    template<class TDistributionModel>
    class CachedSampler : public SimpleSampler<TDistributionModel> {
        friend void testViz();
        template <class TDistributionModel> friend class CacheViz;
        template <class TDistributionModel> friend class RecordViz;
        friend struct State;
        template <class TDistributionModel> friend class CachedSamplerVizAPI;

        ImportanceCache<TDistributionModel> cache;

        Config config;

        bool isCreateNewRecords;
    public:        

        CachedSampler() {            
        }

        virtual bool setCreateNewRecords( bool val ) { 
            isCreateNewRecords = val; 
            if ( !isCreateNewRecords ) {
                ILog( EWarn, "Creating new records in the cache was forbidden." );
            }

            return true;
        }

        virtual VizAPI * getVizApi() { return new CachedSamplerVizAPI<TDistributionModel>( this ); }

        const IStack<typename ImportanceCache<TDistributionModel>::Record*>& getRecords() const {
            return cache.getRecords();
        }

#ifdef LIBIMP_GATHER_PARTICLES
        void loadGatherPoints() {
            SimpleSampler<TDistributionModel>::loadGatherPoints();
            ResultBuffer buffer;
            for ( auto gp = m_gatherPoints.begin(); gp < m_gatherPoints.end(); ++gp ) {                
                Distribution * dist = getDistribution( *gp, buffer );
                if ( dist ) {
                    dist->release();
                }
            }
        }
#endif

        virtual void initImpl(InputIterator& samples, const Config& _config, const Camera* camera, Stats* stats) {
            isCreateNewRecords = true;
            this->config = _config;
            SimpleSampler<TDistributionModel>::initImpl(samples, this->config, camera, stats);
            const BoundingBox3 bbox = BoundingBox3(this->config.sceneBboxMin, this->config.sceneBboxMax);
            cache.init(this->config, stats, bbox);
        }


        float getInterpolationMetric(const Hit& hit, const InterpolationMetricType type) const {
            if(type == METRIC_KL_DIVERGENCE) {
                ResultBuffer tmpBuffer;
                float res;
                cache.getInterpolated(hit, 10.f, config->cache.maxRadius, tmpBuffer, res, false);
                return res/6.f;
            } else if(type == METRIC_IRRADIANCE) {
                ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> buffer;

                const int found = this->nnquery(hit, buffer);
                if(found == 0) {
                    return 0.f;
                }
                const float radius = std::sqrt(buffer[0].weight);
                float result = 0.f;
                for(int i = 0; i < found; ++i) {
                    const Particle& actual = this->particles[buffer[i].index];
                    if(dot(actual.incidentDir, hit.normal) < 0.f) {
                        const float distance = (actual.position-hit.position).size();
                        result += actual.weight*actual.getWeight(distance, radius);
                    }
                }
                result /= 2*PI*radius*radius;

                return result;
            } else if(type == METRIC_D_IRRADIANCE) {
                ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> buffer;

                const int found = this->nnquery(hit, buffer);
                if(found == 0) {
                    return 0.f;
                }
                const float radius = std::sqrt(buffer[0].weight);
                float irradiance = 0.f;
                float dX = 0.f;
                float dY = 0.f;
                Vector3 xDir, yDir;
                coordinateSystem(hit.normal, xDir, yDir);
                IMPORTANCE_ASSERT(xDir.isNormalized() && yDir.isNormalized());
                for(int i = 0; i < found; ++i) {
                    const Particle& actual = this->particles[buffer[i].index];
                    if(dot(actual.incidentDir, hit.normal) < 0.f) {
                        const float distance = (actual.position-hit.position).size();
                        irradiance += actual.weight*actual.getWeight(distance, radius);
                        const Vector3 normDir = (actual.position-hit.position).getNormalized();
                        dX += actual.weight*dot(normDir, xDir)/radius;
                        dY += actual.weight*dot(normDir, yDir)/radius;
                    }
                }
                irradiance /= 2*PI*radius*radius;
                dX /= 2*PI*radius*radius;
                dY /= 2*PI*radius*radius;

                return 3*Vector3(dX, dY, 0).size()/irradiance;
            } else if(type == METRIC_TEST) {
                return abs(hit.normal.x)/2 + abs(hit.normal.y)/2;            
            } else {
                IMPORTANCE_ASSERT(false);
                return 0.f;
            }
        }

        // vraci NULL, kdyz nebylo vzorkovani uspesny
        virtual Distribution* getDistributionImpl( const Hit& hit, IResultBuffer& buffer ) {
            const Float pixel2world = config.pixel2WorldRatio(hit, camera);
            float dummy;
            InterpolationResult<TDistributionModel>* tmpResult = cache.getInterpolated(hit, pixel2world, config.cache.maxRadius, buffer, dummy);
            if(tmpResult) {
                return tmpResult;
            }

            // fallthrough for unsuccessful interpolation
            if ( !isCreateNewRecords ) {
                return NULL;
            }

            ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> particles;   
            const int found = this->nnquery(hit, particles);
            IMPORTANCE_ASSERT(found < 2 || particles[0].distSqr >= particles[1].distSqr);           
            if ( found == 0 ) {
                return NULL;
            }

            auto * rec = cache.createRecord( hit );          
            CacheStats & cstats = rec->distr.m_cacheStats;
            cstats.particleSearchCount  = found;
            cstats.sqrSearchRadius      = particles[ 0 ].distSqr;

            if(TDistributionModel::getDistribution(distFactory, hit, particles.ptr(), found, rec->distr)) {
                /* If there are no particles left after the filter is applied 
                   (e.g. particles on the other side of the surface are removed)
                   do not add any into the cache. */
                if ( rec->distr.m_cacheStats.particleSearchCount == 0 ) {
                    delete rec;
                    return NULL;
                }
                const TDistributionModel* inserted = 
                            cache.add(rec, pixel2world, config.cache.maxRadius);                
                return new((void*)&buffer) CrossfadeInterpolationResult<TDistributionModel>(inserted);
            } else {
                // factory of any distribution representation should always return something
                // at least it should degrade to cosine or uniform sampling over hemisphere
                // sampler itself should decide what will do with the distribution here
                // it can return NULL when there is not enough particles but in that case it must 
                // release the distribution in order to avoid memory leaking                
                IMPORTANCE_ASSERT(false);
                delete rec;
                return NULL;
            }
        }

        virtual void refreshSamples( InputIterator & samples, bool refineCache = true ) {            
            IMPORTANCE_ASSERT( this->config.particles.knn <= MAX_KNN_PARTICLES );            
            tree = PointKdTree<Kd3PositionTraits>();
            particles.clear();

            fill( particles, samples );
            tree.build( particles.begin(), particles.end(), PARTICLE_TREE_LEAF_SIZE );

#ifdef LIBIMP_GATHER_PARTICLES
            gatherParticles();
#endif
			
            int counter = 0;
            auto & recs = cache.getRecords();
            
            #pragma omp parallel for
            for ( int i = 0; i < (int) recs.size(); ++i ) {
                CacheRecord<TDistributionModel> * rec = recs[ i ];
                counter++;
                ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> searchResults;
                Hit hit;
                hit.normal      = rec->normal;
                hit.position    = rec->position;
                const int nFound = nnquery( hit, searchResults );
                CacheStats & cstats         = rec->distr.m_cacheStats;
                cstats.particleSearchCount  = nFound;
                cstats.sqrSearchRadius      = searchResults[ 0 ].distSqr;
                TDistributionModel::getDitributionUpdated( distFactory, hit, searchResults.ptr(), nFound, rec->distr);
				if ( refineCache /*&& nFound > 0*/ ) {
					const Float pixel2world = config.pixel2WorldRatio(hit, camera);
					cache.updateRecordRadius( rec, pixel2world, 
                                              config.cache.maxRadius );					
				}
            }
			if ( refineCache ) {
				cache.rebuildTree();
                ILog( EInfo, "Cache was refined -i.e. we change records radii." );
			}
            ILog( EInfo, "Cache updated - approx. %d records were updated", (int) counter );
        }

        virtual bool isCacheUsed() const {
            return true;
        }
    };
}

#pragma warning(pop)

#endif