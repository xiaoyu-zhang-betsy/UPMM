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


#ifndef ___SIMPLE_SAMPLER___
#define ___SIMPLE_SAMPLER___

#include "../LibImportanceTypes.h"
#include "../sampler/Sampler.h"
#include "../shared/PointKdTree.h"
#include "../viz/vizapi.h"
#include "../shared/particlefilter.h"

#pragma warning(push)
#pragma warning(disable:4265)   //class has virtual functions, but destructor is not virtual
#pragma warning(disable:4239)

namespace Importance {
    
    template<class TDistributionModel>
    class ImmediateDistribution : public Distribution {
        ImmediateDistribution& operator=(const ImmediateDistribution&);
    public:

        typename TDistributionModel model;

        IMPORTANCE_INLINE ImmediateDistribution() {
            IMPORTANCE_ASSERT((size_t(this) % __alignof(ImmediateDistribution)) == 0);
        }

        IMPORTANCE_INLINE ImmediateDistribution(const TDistributionModel& fittingResult) : model(fittingResult) {
            IMPORTANCE_ASSERT((size_t(this) % __alignof(ImmediateDistribution)) == 0);
        }

        virtual void pdfs(const Vector3 * directions, float * pdfs, const int count) const {
            for(int i = 0; i < count; i++) {
                pdfs[i] = float(model.pdf(Vector3(directions[i])));
            }
        }

        virtual void sampleDirections(const Vector2* random, Vector3 * output, float* pdfs, const int count) {
            for(int i = 0; i < count; i++) {
                output[i] = model.sampleDirection(Vector2(random[i])).singlePrecision();
                pdfs[i] = float(model.pdf(Vector3(output[i])));
            }
        }   

		virtual Float gatherAreaPdf(Vector3 wo, Float radius, std::vector<Vector2> &componentCDFs, std::vector<Vector2> &componentBounds){
			IMPORTANCE_ASSERT(false);
			return 1.f;
		}


        virtual void release() {
            model.release();
        }

        virtual const EmInfo * getInfo() const {
            return model.getInfo();
        }

        const void * getModel() const {
            return &model;
        }

        virtual std::string toString() const {
            return model.toString();
        }
    };


    //////////////////////////////////////////////////////////////////////////

    /* Use this in a point kdtree over particles to query for nn particles if you want to ignore particles 
       which come from the opposite side of a surface or those which hit surface "over the corner" */
    class Kd3PositionParticlePredicate {    
    public:
        Kd3PositionParticlePredicate( const IStack<Particle> & particles, const Vector3 & hitNormal ) : m_particles( particles ), m_hitNormal( hitNormal ) {}
        bool operator()( unsigned int index ) const {
            IMPORTANCE_ASSERT( index >= 0 && index < m_particles.size() );
            const Particle & p = m_particles[ index ];
            return ParticleFilterPredicate()( -p.incidentDir, p.normal, m_hitNormal, false );
        }
    private:
        void operator=(const Kd3PositionParticlePredicate & ) {}
        const IStack<Particle> & m_particles;
        const Vector3 m_hitNormal;
    };


    //////////////////////////////////////////////////////////////////////////

        
    template<class TDistributionModel>
    class SimpleSampler : public Sampler {
        friend void testViz();
        template <class TDistributionModel> friend class CacheViz;
        template <class TDistributionModel> friend class RecordViz;
        template <class TDistributionModel> friend class SamplersFacade;
        friend struct State;
    protected:

        IStack<Particle> particles;
        PointKdTree<Kd3PositionTraits> tree;        
        Config config;
        const Camera* camera;        
        Stats* stats;
#ifdef LIBIMP_STATS
        SampleStats m_particleStats;
#endif

    public:
        typedef TDistributionModel DistributionType;
    public:
#ifdef LIBIMP_GATHER_PARTICLES
        std::vector<Hit> m_gatherPoints;
        std::vector<std::vector<Importance::Particle>> m_gatheredParticles;

        void saveGatherPoints( const std::vector<Hit> & hits ) {
            std::ofstream of( "gatherpoints.dump", ios_base::out | ios_base::trunc | ios_base::binary );
            size_t n = hits.size();
            of.write( (char*) &n, sizeof(size_t) );
            for ( auto hit = hits.begin(); hit < hits.end(); ++hit ) {
                of.write( (char*) &hit->position.x, sizeof( float ) );
                of.write( (char*) &hit->position.y, sizeof( float ) );
                of.write( (char*) &hit->position.z, sizeof( float ) );
                of.write( (char*) &hit->normal.x, sizeof( float ) );
                of.write( (char*) &hit->normal.y, sizeof( float ) );
                of.write( (char*) &hit->normal.z, sizeof( float ) );                   
            }
            of.close();
        }

        void loadGatherPoints() {
            std::ifstream ifile( "gatherpoints.dump", ios_base::binary );
            if ( ifile.fail() ) {
                ILog( EWarn, "gatherpoints.dump file can not be read." );
                return;
            }
            size_t n = 0;
            ifile.read( (char*) &n, sizeof(size_t) );
            m_gatherPoints.resize( n );
            size_t i = 0;
            while( ifile && i < n ) {
                Hit & hit = m_gatherPoints[ i ];
                ifile.read( (char*) &hit.position.x, sizeof( float ) );
                ifile.read( (char*) &hit.position.y, sizeof( float ) );
                ifile.read( (char*) &hit.position.z, sizeof( float ) );
                ifile.read( (char*) &hit.normal.x, sizeof( float ) );
                ifile.read( (char*) &hit.normal.y, sizeof( float ) );
                ifile.read( (char*) &hit.normal.z, sizeof( float ) );    
                i++;
            }
            IMPORTANCE_ASSERT( i == n );
            ifile.close();
        }

        void gatherParticles() {
            m_gatheredParticles.resize(m_gatherPoints.size());
            int i = 0;
            for ( auto gp = m_gatherPoints.begin(); gp < m_gatherPoints.end(); ++gp, ++i ) {
                ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> searchResults;
                const int found = nnquery( *gp, searchResults );

                /// We just need to fill in the oUsedParticles stack 
                class Wrapper {
                private:
                    std::vector<Particle> & m_data;
                public:
                    Wrapper( std::vector<Particle> & data ) : m_data( data ) {}
                    void add( const Particle & p, const Vector3 & normal, Float ) {                              
                        Particle tmp = p;
                        // get back to incident direction
                        //tmp.incidentDir = -tmp.incidentDir;
                        m_data.push_back( tmp );
                    }
                } wrapper( m_gatheredParticles[ i ] );      

                DataFilter<Wrapper> filter( &particles, config );
                CacheStats dummy;
                filter.prepareData( searchResults.ptr(), *gp, found, dummy, wrapper );
            }
        }
#endif

        virtual void initImpl(InputIterator& samples, const Config& config, const Camera* camera, Stats* stats) {
            IMPORTANCE_ASSERT(config.particles.knn <= MAX_KNN_PARTICLES);
            this->stats = stats;
            this->config = config;
            this->camera = camera;
            fill( particles, samples );
            distFactory->init(particles, &this->config);
            tree.build(particles.begin(), particles.end(), PARTICLE_TREE_LEAF_SIZE);
#ifdef LIBIMP_GATHER_PARTICLES
            gatherParticles();
#endif
        }

        virtual Distribution* getDistributionImpl(const Hit& hit, IResultBuffer& buffer) {            
            ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> searchResults;
            const int found = nnquery( hit, searchResults );

            ImmediateDistribution<TDistributionModel>* result = createImmediateDistribution( buffer );

           
            const bool success = TDistributionModel::getDistribution(distFactory, hit, searchResults.ptr(), found, result->model);
            
            // factory of any distribution representation should always return something
            // at least it should degrade to cosine or uniform sampling over hemisphere
            // sampler itself should decide what will do with the distribution here
            // it can return NULL when there is not enough particles but in that case it must 
            // release the distribution in order to avoid memory leaking
            IMPORTANCE_ASSERT( success );

            return success ? result : NULL;
        }

        virtual std::string toString() const {
#ifdef LIBIMP_STATS
            std::ostringstream ostr;
            ostr << distFactory->toString() << std::endl
                 << "All particle statistics: " << std::endl
                 << m_particleStats.toString() << std::endl;            
            return ostr.str();
#endif
            return distFactory->toString();
        }

        virtual VizAPI * getVizApi() { return new SimpleSamplerVizAPI<TDistributionModel>( this ); }

        virtual bool isCacheUsed() const {
            return false;
        }

    protected:
        ImmediateDistribution<TDistributionModel> * createImmediateDistribution( IResultBuffer & buffer ) const {
            return new ((void*)&buffer) ImmediateDistribution<TDistributionModel>();
        }

        void fill(IStack<Particle> & particles, InputIterator & samples) {    
            while( samples.isValid() ) {
                Particle photon;
                photon.incidentDir  = Vector3(samples.incidentDirection());
                photon.position     = Vector3(samples.position());
                if ( !config.isWeightOverride ) {
                    photon.weight       = samples.weight();
                } else {
                    photon.weight = 1.f;
                }
				IMPORTANCE_ASSERT( photon.weight > 0.0f );
                photon.normal       = Vector3(samples.normal());
                photon.distance     = samples.distance();
                particles.push(photon);
                samples.next();                
            }
#ifdef LIBIMP_STATS
            computeSampleStatistics( particles, m_particleStats );
#endif
        }

        // pouzivam i v CachedSampler, nerozkopat prosim ;)
        int nnquery( const Hit & hit, ImStaticArray<KdQueryResult, MAX_KNN_PARTICLES> & result ) const {  
            const Float pixel2world = config.pixel2WorldRatio(hit, camera);
            return tree.rangeQuery( 
                Vector3( hit.position ), 
                config.cache.maxRadius * pixel2world, 
                result.ptr(),                
                config.particles.knn/*, Kd3PositionParticlePredicate( particles, hit.normal )*/ );  
        }
    };    
}

#pragma warning(pop)

#endif