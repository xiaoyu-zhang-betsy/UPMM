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

namespace Importance {

    typedef struct stats {
        size_t runningTimeTotal;
        size_t runningTimeEM;
        size_t nConstructionFail;
        size_t nEMRuns;
        size_t nEMIterations;
        __int64 nPhotonsFound;
        __int64 nNNQuery;
        __int64 nDisposedComponents;
        __int64 avgComponents;
        Float avgMaxLikelihood;
        struct stats() {
            nConstructionFail       = 0;
            nEMRuns                 = 0;
            nEMIterations           = 0;
            nPhotonsFound	        = 0;
            nNNQuery                = 0;
            nDisposedComponents		= 0;
            runningTimeTotal		= 0;
            runningTimeEM			= 0;
            avgComponents			= 0;
            avgMaxLikelihood		= 0;
        }
    } EMStats;
}