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
#include "../shared/Vector3.h"

namespace Importance {

class Sphere {
public:
    Sphere() {
        m_radius = -1.0f;
        m_center = Vector3( 0.f );
    }

	Sphere( const Vector3 & center, Float radius, bool flipNormals = false ) {
        IMPORTANCE_ASSERT( radius > 0.f );
		/// Are the sphere normals pointing inwards? default: no
		m_flipNormals   = flipNormals;
        m_radius        = radius;
		m_center        = center;				
    }

    const Vector3 & getCenter() const {
        return m_center;
    }

	bool rayIntersect( const Vector3 & origin, const Vector3 & direction, Float mint, Float maxt, Float & t, Vector3 & outNormal, Vector3 & itsPoint ) const {
		Vector3d o = Vector3d(origin) - Vector3d(m_center);
		Vector3d d(direction);

		double A = d.square();
		double B = 2 * dot(o, d);
		double C = o.square() - m_radius*m_radius;

		double nearT, farT;
		if (!solveQuadraticDouble(A, B, C, nearT, farT))
			return false;

		if (nearT > maxt || farT < mint)
			return false;
		if (nearT < mint) {
			if (farT > maxt)
				return false;
			t = (Float) farT;
		} else {
			t = (Float) nearT;
		}

        itsPoint = origin + t * direction;
        outNormal = (itsPoint - m_center).getNormalized();        
        if ( m_flipNormals ) {
            outNormal = -1.f * outNormal;
        }

		return true;
	}

    Float getRadius() const {
        return m_radius;
    }

    static bool solveQuadraticDouble(double a, double b, double c, double &x0, double &x1) {
	    /* Linear case */
	    if (a == 0) {
		    if (b != 0) {
			    x0 = x1 = -c / b;
			    return true;
		    }
		    return false;
	    }

	    double discrim = b*b - 4.0f*a*c;

	    /* Leave if there is no solution */
	    if (discrim < 0)
		    return false;

	    double temp, sqrtDiscrim = std::sqrt(discrim);

	    /* Numerically stable version of (-b (+/-) sqrtDiscrim) / (2 * a)
	     *
	     * Based on the observation that one solution is always
	     * accurate while the other is not. Finds the solution of
	     * greater magnitude which does not suffer from loss of
	     * precision and then uses the identity x1 * x2 = c / a
	     */
	    if (b < 0)
		    temp = -0.5f * (b - sqrtDiscrim);
	    else
		    temp = -0.5f * (b + sqrtDiscrim);

	    x0 = temp / a;
	    x1 = c / temp;

	    /* Return the results so that x0 < x1 */
	    if (x0 > x1)
		    std::swap(x0, x1);

	    return true;
    }

private:
	Vector3 m_center;
	Float m_radius;
	bool m_flipNormals;
};

}