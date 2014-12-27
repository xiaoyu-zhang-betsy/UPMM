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
#include "../shared/Config.h"

namespace Importance {

struct Float1Traits {
    typedef float TPixel;

    IMPORTANCE_INLINE bool isDifferent(const TPixel x, const TPixel y) const {
        return x != y;
    }
    static IMPORTANCE_INLINE float toScalar(const TPixel x) {
        return x;
    }
    IMPORTANCE_INLINE TPixel zeroPixel() const {
        return 0.f;
    }

    static IMPORTANCE_INLINE TPixel ONE() {
        return 1.f;
    }
};


template<class Traits>
class ImportanceSampler1d  {
    typedef typename Traits::TPixel TPixel;

    Traits traits;

    struct Region {
        float start;
        float size;
        TPixel pdf;

        IMPORTANCE_INLINE float regionProbablity() const {
            return size*Traits::toScalar(pdf);
        }
    };

    IStack<Region> regions;

    IStack<float> cdf;

    IStack<int> lookupTable;

    float lookupStepInv;

    static const int LOOKUP_MULT = 4;

public:

    // vraci avgLuminance
    template<typename Iter>
    TPixel construct(const Iter begin, const Iter end, const Traits& traits) {
        this->traits = traits;
        TPixel avgLuminance = traits.zeroPixel();

        TPixel last = *begin;
        const int count = int(end-begin);
        
        float beginCoord = 0.f;
        for(Iter actual = begin + 1; actual != end+1; ++actual) {
            if(actual == end || traits.isDifferent(*actual, last)) {
                const float coord = int(actual-begin) / float(count);
                IMPORTANCE_ASSERT(coord <= 1.f);
                Region tmp;
                tmp.start = beginCoord;
                tmp.size = coord-beginCoord;
                tmp.pdf = last;
                regions.push(tmp);
                beginCoord = coord;
                avgLuminance += tmp.size*tmp.pdf;
                if(actual != end) {
                    last = *actual;
                }      
            }
        }

        if(regions.size() == 0) {
            Region region;
            region.pdf = Traits::ONE();
            region.start = 0.f;
            region.size = 1.f;
            this->regions.push(region);
        }

        this->cdf.resize(regions.size());
        this->cdf[0] = regions[0].regionProbablity();

        for(int i = 1; i < regions.size(); i++) {
            this->cdf[i] = regions[i].regionProbablity() + cdf[i-1];
        }

        lookupTable.resize(LOOKUP_MULT*cdf.size()+1);
        lookupStepInv = float(lookupTable.size()-1)/cdf[cdf.size()-1];

        int index = 0;
        for(int i = 0; i < lookupTable.size(); i++) {
            const float value = i/lookupStepInv;

            while(index < cdf.size() && cdf[index] < value) {
                index++;
            }
            lookupTable[i] = index;
        }
        lookupTable.end()[-1] = int(cdf.size());
        return avgLuminance;
    }

    IMPORTANCE_INLINE float sample(float random, TPixel& pdf) const {
        const float dummy = random * cdf[cdf.size()-1];

        const int lookupIdx = int(dummy*lookupStepInv);

        //TODO!!!!!!!! - tuhle magii vyresit
        const int from = std::max(0, lookupTable[lookupIdx]-1);
         const int to   = lookupTable[std::min( (size_t) lookupIdx+1, lookupTable.size() - 1 )];

        IMPORTANCE_ASSERT((from == 0 || cdf[from] <= dummy) && (to == cdf.size() || cdf[to] >= dummy));

        const auto result = std::upper_bound(cdf.cbegin()+from, cdf.cbegin()+to, dummy);
        IMPORTANCE_ASSERT(result != cdf.cend());
        const int index = int(result-cdf.cbegin());
        const float previous = (result == cdf.cbegin()) ? 0.f : result[-1];
        random = std::min(1-FLT_EPSILON, (dummy - previous)/(*result - previous));
        
        const Region& region = regions[index];
        const float coords = region.start + random*region.size;

        pdf = region.pdf;
        return coords;
    }

    IMPORTANCE_INLINE int regionsCount() const {
        return this->regions.size();
    }
};

}