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

class VizTriangle {
public:
    Vector3f v0, v1, v2;
    Vector3f color;
    Vector3f n0, n1, n2;

    IMPORTANCE_INLINE float intersect(const Vector3f origin, const Vector3f direction, Vector3f& outNormal, Vector3f& outColor) const {
        const Vector3f pointOrigin = v0 - origin;
        const Vector3f areaNormal = cross(v2-v0, v1-v0);
        const float normalDirDot = dot(areaNormal, direction);
        const float invVolume = 1.f/normalDirDot;
        const float rayParameter = dot(pointOrigin, areaNormal)*invVolume;
        if(rayParameter < 0) {
            return INFINITY;
        }
        const Vector3f _cross = cross(pointOrigin, direction);
        const float volumeY = -dot(_cross, v2-v0);
        const float v = volumeY*invVolume;
        if(v < 0) {
            return INFINITY;
        }

        const float volumeZ = dot(_cross, v1-v0);
        const float w = volumeZ*invVolume;
        const float u = 1.f - v - w;
        if(u < 0 || w < 0) {
            return INFINITY;
        }
        outColor = this->color;
        outNormal = (n0*u + n1*v + n2*w).getNormalized();
        if(dot(outNormal, direction) > 0.f) {
            outNormal *= -1;
        }
        return rayParameter;
    }
};


class TriangleIterator {
public:
    // returns true if there was another triangle
    virtual bool next(Vector3f& outVert0, Vector3f& outVert1, Vector3f& outVert2, 
        Vector3f& outNormal0, Vector3f& outNormal1, Vector3f& outNormal2, Vector3f& outColor) = 0;
};

}