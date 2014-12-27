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

#include "..\shared\Utils.h"
#include  "..\shared\Particle.h"
#include "..\shared\Stack.h"
#include "..\shared\Bitmap.h"

namespace Importance {
    class EnviroMap {
    public:
        virtual ~EnviroMap() { }
        virtual float _dirPdf(const Vector3 dirFromEnviro) const = 0;
        virtual Vector3 getSceneCenter() const = 0;        
        virtual float getDiscRadius() const = 0;
    };

    namespace EnviroSamplerMapping {
        inline Vector3 getDir(const Vector2 uniformSquare) {
            IMPORTANCE_ASSERT(uniformSquare.x >= 0.f && uniformSquare.x <= 1.f && uniformSquare.y >= 0.f && uniformSquare.y <= 1.f);
            return uniformSphere(uniformSquare.x, uniformSquare.y);
        }

        inline Vector2 getDirInverse(const Vector3 in) {
            Vector2 result;
            uniformSphereInverse(in, result.x, result.y);
            IMPORTANCE_ASSERT( result.isValid() );
            return result;
        }
    }

    class DebuggingBundle {
    public:
        DebuggingBundle() {
            m_all[ 0 ] = &m_counts;
            m_all[ 1 ] = &m_hist;
            m_all[ 2 ] = &m_pdf;
            m_all[ 3 ] = &m_totalCounts;

            m_names[ 0 ] = "counts";
            m_names[ 1 ] = "histogram";
            m_names[ 2 ] = "pdf";
            m_names[ 3 ] = "totalCounts";
        }

        IMPORTANCE_INLINE void updateHistogramAndCounts( const IStack<Particle> & particles ) {
            const Vector2 dims( (Float) m_hist.getWidth(), (Float) m_hist.getHeight() );
            m_counts.fill( 0.f );
            for ( auto it = particles.cbegin(); it < particles.cend(); ++it ) {
                const Vector2 point = EnviroSamplerMapping::getDirInverse( it->position ) * dims;
                m_counts( std::min( (int) point.x, m_counts.getWidth()-1 ), 
                    std::min( (int) point.y, m_counts.getHeight()-1 ) ) += 1.f;
                m_totalCounts( std::min( (int) point.x, m_totalCounts.getWidth()-1 ), 
                               std::min( (int) point.y, m_totalCounts.getHeight()-1 ) ) += 1.f;
                m_hist( std::min( (int) point.x, m_hist.getWidth()-1 ),
                        std::min( (int) point.y, m_hist.getHeight()-1 ) ) += it->weight;
            }
        }

        IMPORTANCE_INLINE void resize( int x, int y ) {
            for ( auto it = m_all.begin(); it < m_all.end(); ++it ) {
                (*it)->resize( x, y );
            }
        }

        IMPORTANCE_INLINE Float sumHistogram() {
            Float sum = 0.f;
            for ( int x = 0; x < m_hist.getWidth(); ++x ) {
                for ( int y = 0; y < m_hist.getHeight(); ++y ) {
                    sum += m_hist( x, y );
                }
            }
            return sum;
        }

        IMPORTANCE_INLINE void normalize() {
            Float sum = sumHistogram();
            IMPORTANCE_ASSERT( sum > 0.f );
            Float invSum = 1.f / sum;
            for ( int x = 0; x < m_hist.getWidth(); ++x ) {
                for ( int y = 0; y < m_hist.getHeight(); ++y ) {
                    m_hist( x, y ) *= invSum;
                }
            }
        }
    public:
        ImBitmap<Float> m_totalCounts;
        ImBitmap<Float> m_counts;
        ImBitmap<Float> m_hist;
        ImBitmap<Float> m_pdf;                

        ImStaticArray<ImBitmap<Float>*, 4> m_all;
        ImStaticArray<std::string, 4> m_names;
    };
}