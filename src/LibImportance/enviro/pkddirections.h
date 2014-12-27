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
#include "../shared/PointKdTree.h"
#include "progressivedensity.h"
#include "ienviromap.h"

namespace Importance {

    /* Used for clamping. Zero means no clamping. */
    const float MIN_PARTICLE_WEIGHT             = 0.f;//1e-5f;


    class PKDDirections : public IDensityEstimator {
    public:
        PKDDirections( const IStack<Particle> & particles, const PointKdTree<Kd3PositionTraits> & tree ) 
              : m_densityEstimator( NULL ), 
                m_totalWeightInv( 0.f ),
                m_particles( particles ),
                m_tree( tree ), 
                m_rasterHeight( 0 ) {}


        virtual ~PKDDirections() {
            if ( m_densityEstimator ) {
                delete m_densityEstimator;
            }
        }


        IMPORTANCE_INLINE void init( const EnviroMap * enviro, 
            ImBitmap<Float> & raster, DebuggingBundle & debug, const Config & /*cfg*/ ) {            

            int w = raster.getWidth(),
                h = raster.getHeight();            
            m_rasterHeight = h;

            const Vector2 dims((Float)w, (Float)h);
            m_pixels.reserve( w * h );
            for(int x = 0; x < w; ++x) {
                for(int y = 0; y < h; ++y) {                                     
                    PDensityEst::Pixel pixel;
                    pixel.value     = 0.f;
                    pixel.position  = EnviroSamplerMapping::getDir(Vector2(x+0.5f, y+0.5f)/dims);
                    m_pixels.push( pixel );
                }
            }

            IMPORTANCE_ASSERT( m_densityEstimator == NULL );
            m_densityEstimator = new PDensityEst( 0.7f, *this );            
            densityEstimateDirections( enviro, true, raster, debug );
        }


        IMPORTANCE_INLINE void refresh( const EnviroMap * enviro, 
            ImBitmap<Float> & raster, DebuggingBundle & debug ) {
                m_densityEstimator->nextIteration();
                densityEstimateDirections( enviro, false, raster, debug );
        }


        std::string toString() const {
            std::ostringstream ostr;
            ostr << "Progressive estimates: " << ( m_densityEstimator ? m_densityEstimator->getCountUpdates() : 0 ) << std::endl;
            return ostr.str();
        }

        //////////////////////////////////////////////////////////////////////////
        // IDensityEstimator interface
        //////////////////////////////////////////////////////////////////////////


        virtual Float computeDensity( const ImStaticArray<KdQueryResult, ENVIRO_DIR_SEARCH_COUNT> & nearest, int found, Float sqrRadius ) const {
            if ( found == 0 ) {
                return 0.f;
            }
            // size of a spherical cap
            Float area               = PI * sqrRadius;
            Float volumeWeight       = averagedWeight( m_particles, nearest, found, sqrRadius );          
            Float importanceDensity  = std::min( volumeWeight / area, std::numeric_limits<Float>::max() );
            importanceDensity       *= m_totalWeightInv;            
            IMPORTANCE_ASSERT( isReal(importanceDensity) && importanceDensity >= 0.f );

            return importanceDensity;
        }



        virtual void knnQuery( const Vector3 & position, int minParticles, int & found, Float & sqrSearchRad, 
            ImStaticArray<KdQueryResult, ENVIRO_DIR_SEARCH_COUNT> & nearest ) const {
                IMPORTANCE_ASSERT( ENVIRO_DIR_SEARCH_COUNT > 0 );            
                const Vector3 & dir = position;
                //square of Maximal Euclid. Distance that we allow
                sqrSearchRad = 2-2*cos( PI/20.f );
                found = m_tree.rangeQuery( dir, std::sqrt( sqrSearchRad ), nearest.ptr(), ENVIRO_DIR_SEARCH_COUNT );
                if ( found >= minParticles ) {
                    sqrSearchRad = nearest[ 0 ].distSqr;
                }
        }

        //////////////////////////////////////////////////////////////////////////            

    protected:
        IMPORTANCE_INLINE void densityEstimateDirections( const EnviroMap* enviro, bool initialEstimate,
            ImBitmap<Float> & raster, DebuggingBundle & debug ) {            
            IMPORTANCE_ASSERT( m_densityEstimator != NULL );
                        
            computeTotalWeight( m_particles );

            for(int x = 0; x < raster.getWidth(); ++x) {
                for(int y = 0; y < raster.getHeight(); ++y) {                  
                    PDensityEst::Pixel & pixel = getPixel( x, y );
                    IMPORTANCE_ASSERT( pixel.position == EnviroSamplerMapping::getDir(Vector2(x+0.5f, y+0.5f) 
/Vector2( (Float)raster.getWidth(), (Float)raster.getHeight() ) ) );
                    m_densityEstimator->estimateKnnRefRad( pixel );
                    const float tmp = getImportanceMapValue( pixel, enviro );
                    raster(x, y) = tmp;
#ifdef LIBIMP_STATS
                    debug.m_pdf(x, y) = pixel.value;
#endif
                    IMPORTANCE_ASSERT(isReal(raster(x, y)) && raster(x, y) >= 0.f);
                }
            }

#ifdef LIBIMP_STATS
            debug.updateHistogramAndCounts( m_particles );
#endif
        }



        IMPORTANCE_INLINE Float getImportanceMapValue( const PDensityEst::Pixel & pixel, const EnviroMap * enviro ) const {
            const float enviroPdf = enviro != NULL ? enviro->_dirPdf( pixel.position ) : 1.f;
            return enviroPdf * pixel.value;
        }


        IMPORTANCE_INLINE void computeTotalWeight( const IStack<Particle> & particles ) {
            m_totalWeightInv = 0.f;
            for ( auto part = particles.cbegin(); part < particles.cend(); ++part ) {
                m_totalWeightInv += std::max( part->weight, MIN_PARTICLE_WEIGHT );
            }
            m_totalWeightInv = 1.f / m_totalWeightInv;
        }
        

        IMPORTANCE_INLINE Float averagedWeight( const IStack<Particle> & particles, 
            const ImStaticArray<KdQueryResult, ENVIRO_DIR_SEARCH_COUNT> & nearest, int found, Float sqrBandwidth ) const {
            if ( found == 0 ) {
                return 0.f;
            }

            Float res = 0.f;            
            // We are using Epanechnikov's kernel
            for ( int i = 0; i < found; ++i ) {
                IMPORTANCE_ASSERT( nearest[ i ].distSqr <= sqrBandwidth );
                Float kernelW = std::max( 0.f, 1.f - nearest[ i ].distSqr / sqrBandwidth );
                res += kernelW * std::max( particles[ nearest[ i ].index ].weight, MIN_PARTICLE_WEIGHT );
            }            
            return res * 3.f / 4.f;
        }


        IMPORTANCE_INLINE PDensityEst::Pixel & getPixel( int x, int y ) {
            IMPORTANCE_ASSERT( m_rasterHeight > 0 );
            size_t index = m_rasterHeight * x + y;
            IMPORTANCE_ASSERT( index < m_pixels.size() );
            return m_pixels[ index ];
        }

    private:
        void operator=( const PKDDirections & ) {};
    public: /* debug */
        IStack<PDensityEst::Pixel> m_pixels;        
    private:
        PDensityEst * m_densityEstimator;
        Float m_totalWeightInv;        
        const IStack<Particle> & m_particles;
        const PointKdTree<Kd3PositionTraits> & m_tree;
        int m_rasterHeight;
    };
}