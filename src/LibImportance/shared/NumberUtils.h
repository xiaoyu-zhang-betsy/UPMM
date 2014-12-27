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

#include <algorithm>
#include "Config.h"
#undef min
#undef max

namespace Importance {

    IMPORTANCE_INLINE bool isNaN(const float x) {        
        int classification = ::_fpclass(x);
        return classification == _FPCLASS_QNAN 
            || classification == _FPCLASS_SNAN;
    }

    IMPORTANCE_INLINE bool isINF(const float x) {
        int type = ::_fpclass(x);
        if (type == _FPCLASS_PINF || type == _FPCLASS_NINF)
            return 1;
        return 0;
    }

    IMPORTANCE_INLINE bool isNaN(const double x) {        
        int classification = ::_fpclass(x);
        return classification == _FPCLASS_QNAN 
            || classification == _FPCLASS_SNAN;
    }

    IMPORTANCE_INLINE bool isINF(const double x) {
        int type = ::_fpclass(x);
        if (type == _FPCLASS_PINF || type == _FPCLASS_NINF)
            return 1;
        return 0;
    }


    /// \brief Returns value clamped to be less or equal to maximum and greater or equal to 
    ///        minimum
    template<class T>
    IMPORTANCE_INLINE T clamp(const T value, const T minimum, const T maximum) {
        IMPORTANCE_ASSERT(minimum <= maximum);
        return std::max(minimum, std::min(value, maximum));
    }

    IMPORTANCE_INLINE bool isReal(const Float x) {
        return x < INFINITY && x > -INFINITY && !isNaN(x);
    }
}