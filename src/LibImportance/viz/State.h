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
#include "..\caching\CachedSampler.h"
#include "..\enviro\EnviroSampler.h"
#include "..\shared\Triangle.h"
#include "sphere.h"

namespace Importance {
    namespace VizTools {
        const Float g_fcPalette[][3] = {
            {1.0f, 0.0f, 0.0f},
            {1.0f, 1.0f, 0.0f},
            {0.0f, 1.0f, 0.0f},
            {0.0f, 1.0f, 1.0f},
            {0.0f, 0.0f, 1.0f}
        };

        //inline void getColor(float in, unsigned char& r, unsigned char& g, unsigned char& b) {
        //    in = std::min(in, 1.f);
        //    r = unsigned char(255 * std::max(1-abs(in-1.f)*2, 0.f));
        //    g = unsigned char(255 * std::max(1-abs(in-0.5f)*2, 0.f));
        //    b = unsigned char(255 * std::max(1-abs(in-0.f)*2, 0.f));
        //}

        inline void getColor2( Float weight, Vector3 * color ) {
            IMPORTANCE_ASSERT( weight >= 0.f && weight <= 1.f );  
            weight = 1.f - weight;
            int i = (int) ( weight * 4.f );        
            i = std::min( 3, i );
            IMPORTANCE_ASSERT( i >=0 && i <= 3 );
            Vector3 a( g_fcPalette[ i ][0], g_fcPalette[ i ][1], g_fcPalette[ i ][2] );
            Vector3 b( g_fcPalette[ i + 1 ][0], g_fcPalette[ i + 1 ][1], g_fcPalette[ i + 1 ][2] );
            Float t = 4.f * weight - i;
            *color = a * ( 1 - t ) + b * t;            
        }
    }

    enum MappingType {
        MAPPING_HEMISPHERE,
        MAPPING_SHIRLEY,
    };

    struct InitCamera {
        Vector3f origin;
        Vector3f target;
        Vector3f roll;
        float xScale;

        InitCamera() : xScale( 1.f ) { }

        InitCamera( const Vector3f & origin, const Vector3f & target, const Vector3f & roll ) {
            this->origin    = origin;
            this->target    = target;
            this->roll      = roll;
            this->xScale    = 1.f;
        }

        InitCamera( const Vector3f & origin, const Vector3f & target, const Vector3f & roll, float xScale ) {
            this->origin    = origin;
            this->target    = target;
            this->roll      = roll;
            this->xScale    = xScale;
        }
    };

    struct SceneStatistics {
        SceneStatistics() : particleMinWeight( 0.0 ), particleMaxWeight( 0.0 ), particleAvgWeight( 0.0 ) {}
        double particleMaxWeight, particleMinWeight, particleAvgWeight;
    };

    struct State {
        enum EPrimType {
            PT_TRIANGLE,
            PT_SPHERE
        };                        

        Sphere enviroDummySph;
        
        bool showPhotons;
        int everyNthPhoton;

        bool showRecords;
        int everyNthRecord;

        bool showGeometry;

        bool isGeometryBgColor;
        
        bool isRecordsSameColor;

        bool isPhotonsFalseColor;        

        Hit lastClick;

        bool showRadii;

        bool interpolateOnClick;

        float photonLength;

        IStack<VizTriangle> triangles;

        SceneStatistics sceneStats;

        bool showActivePhotons;

        bool showCameraRatio;

        Config config;

        InterpolationMetricType displayedMetric;       

        const std::vector<Hit> * marks;

        bool isShowMarks;

        bool isDisplayNormals;

        SampleStats smplStats;

        BoundingBox3 bbox;

        bool isDisplayUsedParticles;

        bool showEnvironmentOnly;

        Vector3 bgColor;

        struct Record {
            bool isEnviro;
            int clickedX;
            int clickedY;
            float clickedPdf;

            bool showPhotons;
            bool showGeometry;
            bool relativeColorMapping;     

			int gamma;

            MappingType mapping;
            
            ResultBuffer wharrgarbl;
            /// tady na tohle si udelej nakej spolecnej iface pro zobrazovany interpolation resulty
            //Distribution * clickedDistr;
            VizHitRecord vRec;
            
            Record() {
                reset();
            }
            void reset() {                
                clickedX                = -1;
                clickedPdf              = NAN;
                showGeometry            = false;
                showPhotons             = true;
                mapping                 = MAPPING_SHIRLEY;
                relativeColorMapping    = true;
                isEnviro                = false;
				gamma = 2.2f;
            }
        } recordViz;

        State() {            
            isPhotonsFalseColor = true;
            showPhotons = false;
            showRecords = true;
            everyNthRecord = 1;
            everyNthPhoton = 1;
            showGeometry = false;
            isGeometryBgColor = false;
            isRecordsSameColor = false;
            showRadii = false;
            photonLength = 1.f;
            interpolateOnClick = false;
            showActivePhotons = true;
            showCameraRatio = false;        
            displayedMetric = METRIC_NONE;
            marks           = NULL;
            isShowMarks     = false;
            isDisplayNormals = false;
            isDisplayUsedParticles = false;
            showEnvironmentOnly = false;
            bgColor = Vector3( 0.2f );
        }

        IMPORTANCE_INLINE float pdf(const Vector3 dir) {
            if(recordViz.vRec.isValid) {
                return recordViz.vRec.pDistr->pdf( dir );
            }

            return NAN;
        }

        //const typename ImportanceCache<TDistributionModel>::Record * getSelectedRecord() const {
        //    IMPORTANCE_ASSERT( recordViz.selectedEnviroId == -1 || recordViz.selectedId == -1 );

        //    if ( recordViz.selectedEnviroId != -1 ) {
        //        return m_enviroSampler->m_cache.records[ recordViz.selectedEnviroId ];
        //    }

        //    if ( recordViz.selectedId != -1 ) {
        //        auto * sampler = getCachedSampler();
        //        if ( sampler != NULL ) {
        //            return getCachedSampler()->cache.records[ recordViz.selectedId ];
        //        }
        //    }

        //    return NULL;
        //}

        void print ( std::ostream & out, const std::string & str ) const {
            out << str;
            out << std::endl;
        }

        std::string getRecordDesc() {
            std::stringstream result;
            //if(interpolateOnClick && recordViz.clickedDistr) {
            //    CrossfadeInterpolationResult<TDistributionModel>* tmp = dynamic_cast<CrossfadeInterpolationResult<TDistributionModel>*>(recordViz.clickedDistr);
            //    if(tmp) {
            //        result << "Clicked PDF: " << recordViz.clickedPdf << "\n";
            //        result << "Interpolated from " << tmp->found << " records:\n";
            //        for(int i = 0; i < tmp->found; ++i) {
            //            result << i << ": weight: "<< tmp->getWeight(i)/ tmp->cdf[tmp->found-1] << "\n";
            //        }
            //    }
            //} else if(recordViz.selectedId >= 0 && getCachedSampler() ) {
            //    result << "Individual Record " << recordViz.selectedId << ", clicked PDF: ";
            //    result << recordViz.clickedPdf << "\n";
            //    result << getCachedSampler()->cache.records[recordViz.selectedId]->distr.toString();
            //}

            result << "getRecordDesc() is not implemented yet" << std::endl;

            result << std::endl;
            //result << "Particle max weight: " << sceneStats.particleMaxWeight << std::endl;
            //result << "Particle min weight: " << sceneStats.particleMinWeight << std::endl;
            //result << "Particle avg weight: " << sceneStats.particleAvgWeight << std::endl;

            return result.str();
        }


        Hit selectedRecordHit() const {
            if(interpolateOnClick) {
                return lastClick;
            } 
            
            return recordViz.vRec.toHit();            
        }

        Hit rayTrace(Vector3 origin, const Vector3 dir, Vector3& actualColor, int & primType ) const {
            origin += dir/100.f;
            Hit result;
            result.position = Vector3(NAN);            
            Vector3 color, actualNormal;

            struct Intersection {
                float dist;
                Vector3 normal, color, position;
                EPrimType primType;
                Intersection() : dist( INFINITY ), position( Vector3(NAN) ) {}
            } gIts, sphIts1, sphIts2;

            if ( !showEnvironmentOnly ) {
                for(int i = 0; i < this->triangles.size(); ++i) {                
                    const float t = this->triangles[i].intersect(origin, dir, actualNormal, color);
                    if(t < gIts.dist) {
                        gIts.dist     = t;
                        gIts.color    = color;
                        gIts.normal   = actualNormal;
                        gIts.position = origin + t * dir;
                        gIts.primType = PT_TRIANGLE;
                    }
                }
            }

            if ( this->enviroDummySph.getRadius() > 0.f ) {
                Float t;
                Vector3 itsPoint;
                if ( this->enviroDummySph.rayIntersect( origin, dir, 1e-4f, INFINITY, t, actualNormal, itsPoint ) ) {
                    sphIts1.dist     = t;
                    sphIts1.color    = Vector3( 1.f, 1.f, 1.f );
                    sphIts1.normal   = actualNormal;
                    sphIts1.position = itsPoint;
                    sphIts1.primType = PT_SPHERE;
                    /* We want an intersection from inside of the environment sphere. */
                    if ( this->enviroDummySph.rayIntersect( 
                        itsPoint, dir, 1e-4f, INFINITY, t, actualNormal, itsPoint ) ) {
                            sphIts2.dist     = t;
                            sphIts2.color    = Vector3( 1.f, 1.f, 1.f );
                            sphIts2.normal   = actualNormal;
                            sphIts2.position = itsPoint;
                            sphIts2.primType = PT_SPHERE;
                    }
                }
            }

            if ( gIts.dist == INFINITY && sphIts1.dist == INFINITY ) {
                return result;
            }

            Intersection * its = NULL;
            if ( gIts.dist < INFINITY ) {
                its = &gIts;
            } else if ( sphIts2.dist < INFINITY ) {
                its = &sphIts2;
            } else if ( sphIts1.dist < INFINITY ) {
                its = &sphIts1;
            }

            if ( its ) {
                result.normal   = its->normal;
                result.position = its->position;
                primType        = its->primType;
                actualColor     = its->color;                
            }

            return result;
        }

        IMPORTANCE_INLINE Vector3 getNormal() const {
            if( recordViz.vRec.isValid ) {
                return recordViz.vRec.normal;
            }

            return Vector3(NAN);            
        }

        IMPORTANCE_INLINE Vector3 coordsToDir(const int x, const int y, const int sizeX, const int sizeY, Vector3 normal ) const {            
            if(!normal.isReal()) {
                return Vector3(NAN);
            }
            Vector2 tmp(x/float(sizeX), float(y)/sizeY);
            Vector3 dir;
            if(recordViz.mapping == MAPPING_HEMISPHERE) {
                tmp = tmp - Vector2(0.5f);
                tmp = tmp * 2;
                dir = diskToHemisphere(tmp);
            } else {
                dir = uniformHemisphereShirley(tmp.x, tmp.y);
            }
            return Frame(normal).toWorld(dir);
        }


        IMPORTANCE_INLINE Vector3 coordsToDir(const int x, const int y, const int sizeX, const int sizeY) const {
            return coordsToDir(x, y, sizeX, sizeY, getNormal() );
        }

        IMPORTANCE_INLINE void dirToCoords(Vector3 dir, int& x, int& y, const int sizeX, const int sizeY) const {
            const Vector3 normal = getNormal();
            if(!normal.isReal()) {
                x = y = 0;
                return;
            }
            dir = Frame(normal).toLocal(dir);
            Vector2 tmp;
            if(recordViz.mapping == MAPPING_HEMISPHERE) {
                tmp = hemisphereToDisk(dir);
                tmp = tmp/2;
                tmp = tmp + Vector2(0.5f);
            } else {
                uniformHemisphereShirleyInverse(dir, tmp.x, tmp.y);
            }
            x = tmp.x*sizeX;
            y = tmp.y*sizeY;
        }

    };
}