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

#include "..\LibImportanceTypes.h"

namespace Importance {

    struct StepwiseEMConfig {

        StepwiseEMConfig( Float alpha, const Config & cfg ) {
            this->alpha                 = alpha;
            this->batchSize             = cfg.fitting.batchSize;
            this->delayedUpdate         = cfg.fitting.delayedUpdate;
            this->regularizationType    = cfg.fitting.regularizationType;

            /// stepwise
            this->nMaxEMIter            = cfg.fitting.nMaxEMIter;
        }

        StepwiseEMConfig() {
            alpha               = 0.6f;
            batchSize           = 1;
            delayedUpdate       = 10;
            regularizationType  = EMAP;
            nMaxEMIter          = std::numeric_limits<size_t>::max();
        }

        Float alpha;
        int batchSize;
        unsigned int delayedUpdate;
        ERegularizationType regularizationType;
        size_t nMaxEMIter;
    };
}