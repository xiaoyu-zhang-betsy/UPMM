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


/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <sstream> 
#include <string>
#include "Utils.h"
#include "QuadVector.h"

namespace Importance {

    /**
     * \brief Stores a three-dimensional orthonormal coordinate frame
     *
     * This class is mostly used to quickly convert between different
     * cartesian coordinate systems and to efficiently compute certain
     * quantities (e.g. \ref cosTheta(), \ref tanTheta, ..).
     *
     * \ingroup libcore
     * \ingroup libpython
     */
    struct Frame {
	    Vector3 s, t;
	    Vector3 n;

	    /// Default constructor -- performs no initialization!
	    inline Frame() { }

	    /// Given a normal and tangent vectors, construct a new coordinate frame
	    inline Frame(const Vector3 &s, const Vector3 &t, const Vector3 &n)
	     : s(s), t(t), n(n) {
	    }

	    /// Construct a new coordinate frame from a single vector
	    inline Frame(const Vector3 &n) : n(n) {
		    coordinateSystem(n, s, t);
	    }

        inline Frame getInverse() const {
            Frame result = *this;
            std::swap(result.s.z, result.n.x);
            std::swap(result.s.y, result.t.x);
            std::swap(result.t.z, result.n.y);
            return result;
        }

        inline Frame operator*(const Frame& other) const {
               const Frame inv = other.getInverse();

                return Frame(Vector3(dot(s, inv.s), dot(s, inv.t), dot(s, inv.n)),
                             Vector3(dot(t, inv.s), dot(t, inv.t), dot(t, inv.n)),
                             Vector3(dot(n, inv.s), dot(n, inv.t), dot(n, inv.n)));        
        }

	    /// Convert from world coordinates to local coordinates
	    inline Vector3 toLocal(const Vector3 &v) const {
		    return Vector3(
			    Importance::dot(v, s),
			    Importance::dot(v, t),
			    Importance::dot(v, n)
		    );
	    }
        inline QuadVector3 toLocal(const QuadVector3 &v) const {
            return QuadVector3(
                QuadVector3::dot(v, s),
                QuadVector3::dot(v, t),
                QuadVector3::dot(v, n)
                );
        }

	    /// Convert from local coordinates to world coordinates
	    inline Vector3 toWorld(const Vector3 &v) const {
		    return s * v.x + t * v.y + n * v.z;
	    }

	    /** \brief Assuming that the given direction is in the local coordinate 
	     * system, return the cosine of the angle between the normal and v */
	    inline static Float cosTheta(const Vector3 &v) {
		    return v.z;
	    }

	    /** \brief Assuming that the given direction is in the local coordinate
	     * system, return the sine of the angle between the normal and v */
	    inline static Float sinTheta(const Vector3 &v) {
		    Float temp = sinTheta2(v);
		    if (temp <= 0.0f)
			    return 0.0f;
		    return std::sqrt(temp);
	    }

	    /** \brief Assuming that the given direction is in the local coordinate
	     * system, return the tangent of the angle between the normal and v */
	    inline static Float tanTheta(const Vector3 &v) {
		    Float temp = 1 - v.z*v.z;
		    if (temp <= 0.0f)
			    return 0.0f;
		    return std::sqrt(temp) / v.z;
	    }

	    /** \brief Assuming that the given direction is in the local coordinate
	     * system, return the squared sine of the angle between the normal and v */
	    inline static Float sinTheta2(const Vector3 &v) {
		    return 1.0f - v.z * v.z;
	    }

	    /** \brief Assuming that the given direction is in the local coordinate 
	     * system, return the sine of the phi parameter in spherical coordinates */
	    inline static Float sinPhi(const Vector3 &v) {
		    Float sinTheta = Frame::sinTheta(v);
		    if (sinTheta == 0.0f)
			    return 1.0f;
		    return clamp(v.y / sinTheta, -1.0f, 1.0f);
	    }

	    /** \brief Assuming that the given direction is in the local coordinate 
	     * system, return the cosine of the phi parameter in spherical coordinates */
	    inline static Float cosPhi(const Vector3 &v) {
		    Float sinTheta = Frame::sinTheta(v);
		    if (sinTheta == 0.0f)
			    return 1.0f;
		    return clamp(v.x / sinTheta, -1.0f, 1.0f);
	    }

	    /** \brief Assuming that the given direction is in the local coordinate
	     * system, return the squared sine of the phi parameter in  spherical
	     * coordinates */
	    inline static Float sinPhi2(const Vector3 &v) {
		    return clamp(v.y * v.y / sinTheta2(v), 0.0f, 1.0f);
	    }

	    /** \brief Assuming that the given direction is in the local coordinate
	     * system, return the squared cosine of the phi parameter in  spherical
	     * coordinates */
	    inline static Float cosPhi2(const Vector3 &v) {
		    return clamp(v.x * v.x / sinTheta2(v), 0.0f, 1.0f);
	    }

	    /// Return a string representation of this frame
	    inline std::string toString() const {
		    std::ostringstream oss;
		    oss << "Frame[" << std::endl
			    << "  s = " << s.toString() << "," << std::endl
			    << "  t = " << t.toString() << "," << std::endl
			    << "  n = " << n.toString() << std::endl
			    << "]";
		    return oss.str();
	    }

    };

}
