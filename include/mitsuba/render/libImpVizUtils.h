/*
    This file is part of a demo implementation of an importance sampling technique
    described in the "On-line Learning of Parametric Mixture Models for Light Transport Simulation"
    (SIGGRAPH 2014) paper.
    The implementation is based on Mitsuba, a physically based rendering system.

    Copyright (c) 2014 by Jiri Vorba, Ondrej Karlik, Martin Sik.
    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once 

/// !!! Include this header before Mitsuba's logger.h file (Log macro mustn't be defined otherwise there is a clash with wxWidgets headers)

//include platform first so that all necessary macros are defined (e.g. _WIN32_WINNT, _CRT_SECURE_NO_WARNINGS, etc.)
#include <mitsuba/core/platform.h>
#include "LibImportance/viz/Viz.h"
#include <mitsuba/render/scene.h>
#include "LibImportance/LibImportance.h"
#include <mitsuba/render/libImpUtils.h>
#include <mitsuba/core/fstream.h>

MTS_NAMESPACE_BEGIN

class VisTriangles : public Importance::TriangleIterator {
public:
    VisTriangles( const std::vector<TriMesh*> & meshes ) : m_meshes( meshes ) {        
        m_meshCounter   = 0;
        m_triCounter    = 0;
    }

    virtual bool next(Importance::Vector3f& outVert0, Importance::Vector3f& outVert1, Importance::Vector3f& outVert2, 
        Importance::Vector3f& outNormal0, Importance::Vector3f& outNormal1, Importance::Vector3f& outNormal2, 
        Importance::Vector3f& outColor) {
            if ( m_meshCounter >= m_meshes.size() ) {
                return false;
            }

            size_t triCount = m_meshes[ m_meshCounter ]->getTriangleCount();
            if ( m_triCounter >= triCount ) {
                m_meshCounter ++;
                m_triCounter = 0;
                return next( outVert0, outVert1, outVert2, outNormal0, outNormal1, outNormal2, outColor );
            }

            Triangle & t        = m_meshes[ m_meshCounter ]->getTriangles()[ m_triCounter++ ];
            const Point & p0    = m_meshes[ m_meshCounter ]->getVertexPositions()[ t.idx[ 0 ] ];
            const Point & p1    = m_meshes[ m_meshCounter ]->getVertexPositions()[ t.idx[ 1 ] ];
            const Point & p2    = m_meshes[ m_meshCounter ]->getVertexPositions()[ t.idx[ 2 ] ];
            outVert0 = IMP_VECTOR3( p0 );
            outVert1 = IMP_VECTOR3( p1 );
            outVert2 = IMP_VECTOR3( p2 );

            const Normal * normals  = m_meshes[ m_meshCounter ]->getVertexNormals();
            Normal n;            
            if ( normals != NULL ) {                
                const Normal & n0       = normals[ t.idx[ 0 ] ];
                const Normal & n1       = normals[ t.idx[ 1 ] ];
                const Normal & n2       = normals[ t.idx[ 2 ] ];
                outNormal0 = IMP_VECTOR3( n0 );
                outNormal1 = IMP_VECTOR3( n1 );
                outNormal2 = IMP_VECTOR3( n2 );
                n = n0 + n1 + n2;
                if (!n.isZero())
                    n /= n.length();
            } else {
                Vector side1(p1-p0), side2(p2-p0);
                n = Normal(cross(side1, side2));                
                if (!n.isZero())
                    n /= n.length();
                outNormal0 = IMP_VECTOR3( n );
                outNormal1 = IMP_VECTOR3( n );
                outNormal2 = IMP_VECTOR3( n );
            }

            Color3 * colors = m_meshes[ m_meshCounter ]->getVertexColors();
            if ( colors != NULL ) {
                const Color3 & s0 = colors[ t.idx[ 0 ] ];
                const Color3 & s1 = colors[ t.idx[ 1 ] ];
                const Color3 & s2 = colors[ t.idx[ 2 ] ];
                Color3 c = ( s0 + s1 + s2 ) / 3.f;            
                outColor.x = c[ 0 ]; outColor.y = c[ 1 ]; outColor.z = c[ 2 ];
            } else {
                n = n + Normal( 1.f, 1.f, 1.f );
                n *= 0.5f;
                outColor = IMP_VECTOR3( n );
            }

            return true;
    }
private:
    size_t m_meshCounter;
    size_t m_triCounter;
    const std::vector<TriMesh*> & m_meshes;
};

inline const std::string cacheStatsToString( const Importance::Stats & cacheStats, const std::string & cacheName ) {
    std::ostringstream str;
    str << std::endl << cacheName << " CacheStatistics[ " << std::endl
        << "Attempted queries: " << cacheStats.queries.attempted << std::endl
        << "Rejected queries: " << cacheStats.queries.rejected << std::endl
        << "Cache records: " << cacheStats.cache.records << std::endl                    
        << "Cache records found: " << cacheStats.cache.recordsFound << std::endl                    
        << "]";
    
    return str.str();
}

inline void statistics_and_visual( 
    const Scene *scene, 
    const GuidingConfig & cfg,
    Importance::Sampler * sampler, 
    Importance::IEnviroSampler * enviroSampler, 
    const Importance::Stats * cacheStats,
    const Importance::Stats * enviroStats,    
    bool importanceFitting = false,
    const std::vector<Importance::Hit> * markers = NULL )
{    
    if ( cfg.m_importance.cache.useCache && cacheStats != NULL ) {        
        SLog( EInfo, "%s", cacheStatsToString( *cacheStats, importanceFitting ? "Importance" : "Radiance" ).c_str() );
    }

    if ( importanceFitting && enviroSampler != NULL && enviroStats != NULL ) {
        SLog( EInfo, "%s", cacheStatsToString( *enviroStats, "Importance environment" ).c_str() );
    }

    if ( sampler != NULL ) {
        SLog( EInfo, "%s", sampler->toString().c_str() );
    }

    if ( enviroSampler != NULL ) {
        SLog( EInfo, "%s", enviroSampler->toString().c_str() );
    }

    if ( cfg.m_mitsuba.showVisualization ) {
        bool isClosed = false;
        ref<const Sensor> sensor = scene->getSensor();
        Point origin     = sensor->getWorldTransform()->operator()( 0.f, Point( 0.f, 0.f, 0.f ) );
        Vector up        = sensor->getWorldTransform()->operator()( 0.f, Vector( 0.f, 1.f, 0.f ) );
        Vector dir       = sensor->getWorldTransform()->operator()( 0.f, Vector( 0.f, 0.f, 1.f ) );

        Float r = scene->getBSphere().radius;                
        Ray ray( origin, Vector( dir ), 0.f );
        Intersection its;
        if ( scene->rayIntersect( ray, its ) ) {
            r = std::min( its.t, r );
        }

        Importance::InitCamera visCam( IMP_VECTOR3( origin ), IMP_VECTOR3( origin + r * dir ), IMP_VECTOR3( up ), -1.f );
        VisTriangles triangles( scene->getMeshes() );                
        Importance::Config tmpConfig = cfg.m_importance;

        showVizWindow( 
            importanceFitting ? "Importance" : "Radiance", 
            sampler,
            importanceFitting ? enviroSampler : NULL,            
            tmpConfig,
            markers,
            true, 
            &triangles, 
            &isClosed, 
            &visCam, 
            true );
    }    
}

ref<Bitmap> impBmpToMtsBmp( const Importance::ImBitmap<Float>& iBmp ) {
    int w = iBmp.getWidth(),
        h = iBmp.getHeight();
    ref<Bitmap> resBmp = new Bitmap( Bitmap::ELuminance, Bitmap::EFloat32, Vector2i( w, h ) );
    Float * dataPtr = resBmp->getFloat32Data();
    for ( int j = 0; j < h; ++j ) {
        for ( int i = 0; i < w; ++i ) {
            *(dataPtr++) = iBmp.get( i, j );
        }
    }
    return resBmp;
}

void writeImportanceEnviroMap( const std::string& fname, const Importance::ImBitmap<float>& imBitmap ) {
    ref<Bitmap> bmp = impBmpToMtsBmp( imBitmap );
    ref<FileStream> fs = new FileStream( fname, FileStream::ETruncWrite );
    bmp->write( Bitmap::EOpenEXR, fs );
}

MTS_NAMESPACE_END