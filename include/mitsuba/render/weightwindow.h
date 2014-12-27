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

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

class WeightWindow {
public:
    WeightWindow() : lowerBound( 0.f ), upperBound( 0.f ), lightWeightScale( 1.f ), eyeWeightScale( 1.f ), m_useAdjoint( false ) {}

    WeightWindow(float _lowerBound, float size, const Spectrum & _lightWeightScale, const Spectrum & _eyeWeightScale) 
        : lowerBound(_lowerBound), upperBound(_lowerBound+size), 
        lightWeightScale(_lightWeightScale), eyeWeightScale(_eyeWeightScale),
        m_useAdjoint( false ) {
            SAssert(lowerBound >= 0.f);
            SAssert(upperBound > lowerBound);
            SAssert( !lightWeightScale.isZero() );
            SAssert( !eyeWeightScale.isZero() );            
    }

    inline bool shouldSplit(const Spectrum & color, bool useAdjoint) const {        
        int index = getMaxBand(color);
        return color[index] > upperBound * getWeightScale(useAdjoint)[index];        
    }

    inline bool shouldSplit( const Spectrum & color ) const {
        return shouldSplit( color, m_useAdjoint );
    }

    template <class T> 
    static inline int argMax(const T x, const T y, const T z) {
        return (x > y) ? ((x > z) ? 0 : 2) : ((y > z) ? 1 : 2);
    }

    static inline int getMaxBand(const Spectrum & color) {
        SAssert( SPECTRUM_SAMPLES == 3 );
        return argMax( color[0], color[1], color[2] );
    }

    /// Returns number of new particles after splitting and sets their color
    inline int split( const Spectrum & color, bool useAdjoint ) const {        
        int index = getMaxBand(color);
        SAssert(color[index] > upperBound * getWeightScale(useAdjoint)[index]);
        float b = upperBound * getWeightScale(useAdjoint)[index];
        int count = (int) std::ceil( color[index] / b );        
        SAssert(color[index] / (Float) count >= lowerBound * getWeightScale(useAdjoint)[index]);
        return count;
    }

    inline Float rrWeight( const Spectrum & color, bool useAdjoint ) const {        
        int index = getMaxBand(color);        
        const Float dummy = lowerBound * getWeightScale(useAdjoint)[index];
        if ( color[index] < dummy ) {
            SAssert(dummy > 0.f);
            return color[index] / dummy;
        }
        return 1.f;
    }

    inline int split( const Spectrum & color ) const {
        return split( color, m_useAdjoint );
    }

    inline Float rrWeight( const Spectrum & color ) const {
        return rrWeight( color, m_useAdjoint );
    }

    inline void pathTracing() { m_useAdjoint = false; }
    inline void lightTracing() { m_useAdjoint = true; }

    inline float getLowerBound() const { return lowerBound; }
    inline float getUpperBound() const { return upperBound; }

    inline const Spectrum & getWeightScale( bool useAdjoint ) const {        
        return useAdjoint ? lightWeightScale : eyeWeightScale;
    }
    inline const Spectrum & getWeightScale() const { return getWeightScale( m_useAdjoint );
    }

    std::string toString() const {
        std::stringstream ostr;
        ostr << "lower bound: " << lowerBound << std::endl
             << "upper bound: " << upperBound << std::endl
             << "current mode: " << ( m_useAdjoint ? "light" : "path" ) << " tracing" << std::endl
             << "current weight scale: " << getWeightScale().toString() << std::endl;
        return ostr.str();
    }
private:    
    /// Weight window lower boundary
    float lowerBound;
    /// Weight window upper boundary
    float upperBound;
    /// This should be the total emitted flux in the scene
    Spectrum lightWeightScale;
    /// This should be 1
    Spectrum eyeWeightScale;
    /// light tracing (true) or path tracing mode (false)
    bool m_useAdjoint;
};

/* Russian roulette encapsulation. Returns the survival probability. 
   param power is an initial particle weight.
*/
Float russianRoulette( const Spectrum & throughput, const Spectrum & power, Float eta,                       
                       const WeightWindow & ww, bool useWeightWindow ) {    
    if ( useWeightWindow )  {        
        return ww.rrWeight( throughput * power /* Current particle weight */ );
    } 

    /** Standard russian roulette implementation that is used in Mitsuba. 
	    
        Try to keep path weights equal to one,
		while accounting for the solid angle compression at refractive
		index boundaries. Stop with at least some probability to avoid
		getting stuck (e.g. due to total internal reflection)    
    */
    return std::min( throughput.max() * eta * eta, (Float) 0.95f );    
}


class ParticleStats {
private:
    int nParticles;
    int minLength;
    int maxLength;
    int nTotalLength;
    int nSplits;
    int nKills;
    int nNewBySplit;
    float minWindowWeight;
    float maxWindowWeight;
    float maxWeight;
    float minWeight;

public:
    ParticleStats() {
        nNewBySplit     = nKills = nSplits = 0;
        minWeight       = minWindowWeight = std::numeric_limits<float>::max();
        maxWeight       = maxWindowWeight = -std::numeric_limits<float>::max();
        nTotalLength    = 0;
        minLength       = std::numeric_limits<int>::max();
        maxLength       = 0;
        nParticles      = 0;
    }

    inline void observe( const Spectrum & color, const Spectrum & emittedColor, int depth ) {
        handleWeight( color, emittedColor );
        handleLength( depth );
    }

    /* Updates statistics. 
       param Color is current weight of a particle.
       param emittedColor is a spectral power that the particle had after it's been emitted.
       (it might be even only approximation used in weight window)
       */
    inline void handleWeight(const Spectrum & color, const Spectrum & emittedColor) {
        minWeight = std::min(minWeight, color.max());
        maxWeight = std::max(maxWeight, color.max());

        const int index = WeightWindow::getMaxBand(color);
        const float tmp = color[index] / std::max( 1e-6f, emittedColor[index] );
        minWindowWeight = std::min(minWindowWeight, tmp);
        maxWindowWeight = std::max(maxWindowWeight, tmp);
    }

    inline void handleLength( int depth ) {
        nTotalLength += depth;
        minLength = std::min( minLength, depth );
        maxLength = std::max( maxLength, depth );
        nParticles++;
    }

    inline void incNewBySplit( int n ) { nNewBySplit += n; }

    inline void incSplit() { nSplits++; }

    inline void incKills() { nKills++; }

    inline void operator +=(const ParticleStats & s) {
        nSplits += s.nSplits;
        nKills += s.nKills;
        nNewBySplit += s.nNewBySplit;
        minWindowWeight = std::min( minWindowWeight, s.minWindowWeight );
        maxWindowWeight = std::max( maxWindowWeight, s.maxWindowWeight );            
        minWeight = std::min( minWeight, s.minWeight );
        maxWeight = std::max( maxWeight, s.maxWeight );
        minLength = std::min( minLength, s.minLength );
        maxLength = std::max( maxLength, s.maxLength );
        nTotalLength += s.nTotalLength;
        nParticles += s.nParticles;
    }

    std::string toString() const {
        std::stringstream ss;
        ss << "Particle splits: " << nSplits << std::endl
            << "Particle kills: " << nKills << std::endl
            << "Particles after splitting: " << nNewBySplit << std::endl
            << "Measured particles in total: " << nParticles << std::endl
            << "Average path length: " << (nParticles > 0 ? nTotalLength / (Float) nParticles : 0) << std::endl 
            << "Max path length: " << maxLength << std::endl
            << "Min particle weight: " << minWeight << std::endl
            << "Max particle weight: " << maxWeight << std::endl
            << "Min particle window weight: " << minWindowWeight << std::endl
            << "Max particle window weight: " << maxWindowWeight<< std::endl;
        return ss.str();
    }
};

MTS_NAMESPACE_END