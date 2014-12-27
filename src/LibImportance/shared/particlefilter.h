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

#include "Vector3.h"

namespace Importance {
    class ParticleFilterPredicate {
    public:
        /* Omit the photons on the opposite site of a surface and behind the corners. */
        /* particleDir is direction from point of incidence toward the source of particle */
        bool operator()( const Vector3 & particleDir, const Vector3 & photonHitNormal, const Vector3 & hitNormal, bool overSphere ) {
            return overSphere || (
                        /* Photons must not came from direction under the surface */
                        Importance::dot( particleDir, hitNormal ) >= 1e-3f && 
                        /* We consider only photons that hit surface that forms 
                        at least 45 degrees with queried surface */
                        Importance::dot( photonHitNormal, hitNormal ) >= 0.707f 
                        );
        }
    };
}