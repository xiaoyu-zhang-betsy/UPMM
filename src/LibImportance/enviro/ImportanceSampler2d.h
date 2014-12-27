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
#include "../shared/Bitmap.h"
#include "../shared/Vector2.h"

namespace Importance {
    struct FloatTraits {
        typedef float TPixel;
    
        IMPORTANCE_INLINE bool alreadyUsed(const TPixel pixel) const {
            return pixel < -1.f;
        }
        IMPORTANCE_INLINE void setUsed(TPixel& pixel) const {
            pixel -= 1e6f;
        }
        IMPORTANCE_INLINE bool isDifferent(const TPixel x, const TPixel y) const {
            return abs(x-y) > 1e-5f;
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
    class ImportanceSampler2d {
        typedef typename Traits::TPixel TPixel;

        Traits traits;
    
        struct Region {
            Vector2 startCoords;
            Vector2 coordsSize;
            TPixel pdf;

            IMPORTANCE_INLINE float regionProbablity() const {
                return coordsSize.x*coordsSize.y*Traits::toScalar(pdf);
            }
        };

        IStack<Region> regions;

        IStack<float> cdf;

        IStack<int> lookupTable;

        float lookupStepInv;

        static const int LOOKUP_MULT = 4;

    public:

        void setTraits(const Traits& traits) {
            this->traits = traits;
        }

        // vraci avgLuminance
        // znici bitmapu
        float construct(ImBitmap<TPixel>& bitmap) {
            float avgLuminance = 0.f;
            float area = 0.f;

            for(int x = 0; x < bitmap.getWidth(); ++x) {
                for(int y = 0; y < bitmap.getHeight(); ) {

                    const TPixel actualValue = bitmap(x, y);
                    if(traits.alreadyUsed(actualValue)) {
                        y++;
                        continue;
                    }
                    IMPORTANCE_ASSERT(!traits.alreadyUsed(actualValue));
                    const int startY = y;

                    do {
                        y++;
                    } while(y < bitmap.getHeight() && !traits.isDifferent(actualValue, bitmap(x, y)));


                    TPixel sum(traits.zeroPixel());
                    int actualX = x;
                    bool kill = false;
                    while(!kill) {
                        for(int yTemp = startY; yTemp < y; yTemp++) {
                            TPixel& temp = bitmap(actualX, yTemp);
                            IMPORTANCE_ASSERT(!traits.alreadyUsed(temp));
                            sum += temp;
                            traits.setUsed(temp);
                        }
                        actualX++;
                        if(actualX == bitmap.getWidth()) {
                            break;
                        }
                        for(int yTemp = startY; yTemp < y; yTemp++) {                        
                            if(traits.isDifferent(actualValue, bitmap(actualX, yTemp))) {
                                kill = true;
                                break;
                            }
                        }
                    }
                    sum /= float((y-startY)*(actualX - x));
                    //CASSERT(!traits.isDifferent(sum, actualValue));
                    Region region;
                    region.pdf = sum;
                    region.startCoords = Vector2(float(x)      /bitmap.getWidth(), float(startY)/bitmap.getHeight());
                    region.coordsSize  = Vector2(float(actualX)/bitmap.getWidth(), float(y)     /bitmap.getHeight()) - region.startCoords;
                    IMPORTANCE_ASSERT(region.startCoords.x<1.f && region.startCoords.y<1.f && region.regionProbablity() >= 0.f);
                    avgLuminance += region.regionProbablity();

                    if(sum != traits.zeroPixel()) {
                        this->regions.push(region);
                    }
                    area += region.coordsSize.x*region.coordsSize.y;
                }
            }

            if(regions.size() == 0) {
                Region region;
                region.pdf = Traits::ONE();
                region.startCoords = Vector2(0.f, 0.f);
                region.coordsSize  = Vector2(1.f, 1.f);
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

        IMPORTANCE_INLINE Vector2 sample(Vector2 random, TPixel& pdf) const {
            const float dummy = random.x * cdf[cdf.size()-1];

            const int lookupIdx = int(dummy*lookupStepInv);

            //TODO!!!!!!!! - tuhle magii vyresit
            const int from = std::max(0, lookupTable[lookupIdx]-1);
            const int to   = lookupTable[lookupIdx+1];

            IMPORTANCE_ASSERT((from == 0 || cdf[from] <= dummy) && (to == cdf.size() || cdf[to] >= dummy));

            const auto result = std::upper_bound(cdf.begin()+from, cdf.begin()+to, dummy);
        
            const int index = int(result-cdf.begin());
            const float previous = (result == cdf.begin()) ? 0.f : result[-1];
            random.x = std::min(1-FLT_EPSILON, (dummy - previous)/(*result - previous));
            IMPORTANCE_ASSERT(isReal(random.x));

            const Region& region = regions[index];
            const Vector2 coords = region.startCoords + random*region.coordsSize;

            pdf = region.pdf;
            return coords;
        }

        IMPORTANCE_INLINE int regionsCount() const {
            return this->regions.size();
        }
    };
}