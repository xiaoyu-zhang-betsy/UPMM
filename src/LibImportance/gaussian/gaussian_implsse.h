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
/* SERVES FOR SPECIFYING SSE GAUSSIAN TEMPLATES */

#include "../shared/Config.h"
#include "../shared/Sse.h"
#include "gaussianfactorysse.h"

namespace Importance {	
    typedef GaussianMultiMix<MAX_FITTED_LOBES, 4, Float4, SquareHemisphereMapping> GaussianMixtureSSEHemisphere;
    typedef GaussianMultiMix<MAX_FITTED_LOBES, 4, Float4, SquareSphereMapping> GaussianMixtureSSESphere;    
	
    typedef GaussianMixtureFactorySIMD<Float4, 4, GaussianMixtureSSEHemisphere> GaussianMixtureFactorySSEHemisphere;
    typedef GaussianMixtureFactorySIMD<Float4, 4, GaussianMixtureSSEHemisphere> GaussianMixtureFactorySSESphere;

	const Float4 GaussianMixtureSSEHemisphere::EPSILON            = Float4(1e-12f);
	const Float4 GaussianMixtureSSEHemisphere::MIN_DETERMINANT    = Float4(1e-12f);		
	const Float4 GaussianMixtureSSEHemisphere::PRIOR_A            = Float4(1e-2f);
	const Float4 GaussianMixtureSSEHemisphere::PRIOR_B            = Float4(5e-4f);
	const Float4 GaussianMixtureSSEHemisphere::PRIOR_V            = Float4(1e-2f);
	
    const Float4 GaussianMixtureSSESphere::EPSILON            = GaussianMixtureSSEHemisphere::EPSILON;
    const Float4 GaussianMixtureSSESphere::MIN_DETERMINANT    = GaussianMixtureSSEHemisphere::MIN_DETERMINANT;	
    const Float4 GaussianMixtureSSESphere::PRIOR_A            = GaussianMixtureSSEHemisphere::PRIOR_A;
    const Float4 GaussianMixtureSSESphere::PRIOR_B            = GaussianMixtureSSEHemisphere::PRIOR_B;
    const Float4 GaussianMixtureSSESphere::PRIOR_V            = GaussianMixtureSSEHemisphere::PRIOR_V;

    const Float4 VectorTraits< Float4 >::initVal		= Float4( 0.0f );    
};