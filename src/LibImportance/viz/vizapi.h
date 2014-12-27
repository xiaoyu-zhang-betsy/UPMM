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
#include "../caching/CacheTypes.h"

namespace Importance {

    template<class TDistributionModel>
    class SimpleSampler;

    template<class TDistributionModel>
    class CachedSampler;

    template<class TDistributionModel>
    class EnviroSampler;

    struct State;

    struct VizHitRecord {
        VizHitRecord() : isValid( false ), pDistr( NULL ), selectedId( -1 ) {}

        Distribution * pDistr;
        Vector3 normal;
        Vector3 position;
        bool isValid;
        int selectedId;        

        Hit toHit() const {
            Hit hit;
            hit.normal = normal;
            hit.position = position;
            return hit;
        }
    };

    //////////////////////////////////////////////////////////////////////////

    class VizAPI {
    public:
        static void drawRadius( const Vector3 position, const Vector3 normal, const float radius, const Vector3 c = Vector3(1,0,0) );

        virtual VizHitRecord getVizHitRecord( const Hit & hit, IResultBuffer & buffer, const State & state ) = 0;
        virtual void draw( const State & state ) const = 0;
    protected:
        template<typename TDistributionModel>
        static void setRecordColor( const Importance::CacheRecord<TDistributionModel> & rec, const State & state ) {
            if ( state.isRecordsSameColor ) {
                glColor3f(1, 0, 0);
                return;
            }          

            if( rec.distr.isUndersampled() ) {
                glColor3f(0, 1, 1);
            } else if ( rec.originalRadius != rec.validRadius() ) {
                glColor3f(1, 1, 0);
            } else {
                glColor3f(1, 0, 0);
            }
        }
    };

    //////////////////////////////////////////////////////////////////////////

    template<class TDistributionModel>
    class SimpleSamplerVizAPI : public VizAPI {
    public:
        SimpleSamplerVizAPI( SimpleSampler<TDistributionModel> * sampler );

        virtual VizHitRecord getVizHitRecord( const Hit & hit, IResultBuffer & buffer, const State & state );
        virtual void draw( const State & /*state*/ ) const { }
    private:
        SimpleSampler<TDistributionModel> * m_sampler;
    };

    //////////////////////////////////////////////////////////////////////////

    template<class TDistributionModel>
    class CachedSamplerVizAPI : public VizAPI {
    public:
        CachedSamplerVizAPI( CachedSampler<TDistributionModel> * sampler );

        virtual VizHitRecord getVizHitRecord( const Hit & hit, IResultBuffer & buffer, const State & state );
        virtual void draw( const State & state ) const;
    private:
        void paintActiveRecords( const State & state ) const;
        void paintAllNonEnviroRecords( const State & state ) const;
    private:
         CachedSampler<TDistributionModel> * m_sampler;
    };

    //////////////////////////////////////////////////////////////////////////

    template<class TDistributionModel>
    class EnviroSamplerVizAPI : public VizAPI {
    public:
        EnviroSamplerVizAPI( EnviroSampler<TDistributionModel> * sampler );

        virtual VizHitRecord getVizHitRecord( const Hit & hit, IResultBuffer & buffer, const State & state );
        virtual void draw( const State & state ) const;
        static Vector3 recordToPosition( const State & state, const CacheRecord<TDistributionModel> * rec );
    private:
        void paintAllEnviroRecords( const State & state ) const;
    private:
        EnviroSampler<TDistributionModel> * m_sampler;
    };

}