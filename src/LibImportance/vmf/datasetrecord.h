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

namespace Importance {

    /**
        * \brief Generic multi-dimensional bounding box data structure
        *
        * Maintains a component-wise minimum and maximum position and provides
        * various convenience functions to query or change them.
        *
        * \tparam T Underlying point data type (e.g. \c TPoint3<float>)
        * \ingroup libcore
        */
    template <typename T> struct TAABB {
	    typedef T                           PointType;
	    typedef typename T::Scalar          Scalar;
	    typedef typename T::VectorType      VectorType;

	    /**
	        * \brief Create a new invalid bounding box
	        *
	        * Initializes the components of the minimum
	        * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	        * respectively.
	        */
	    inline TAABB() {
		    reset();
	    }

	    /// Create a collapsed AABB from a single point
	    inline TAABB(const PointType &p)
		    : min(p), max(p) { }

	    /// Create a bounding box from two positions
	    inline TAABB(const PointType &min, const PointType &max)
		    : min(min), max(max) {
    #if defined(LIBIMP_DEBUG)
		    for (int i=0; i<PointType::dim; ++i)
			    IMPORTANCE_ASSERT(min[i] <= max[i]);
    #endif
	    }

	    /// Equality test
	    inline bool operator==(const TAABB &aabb) const {
		    return min == aabb.min && max == aabb.max;
	    }

	    /// Inequality test
	    inline bool operator!=(const TAABB &aabb) const {
		    return min != aabb.min || max != aabb.max;
	    }

	    /// Clip to another bounding box
	    inline void clip(const TAABB &aabb) {
		    for (int i=0; i<PointType::dim; ++i) {
			    min[i] = std::max(min[i], aabb.min[i]);
			    max[i] = std::min(max[i], aabb.max[i]);
		    }
	    }

	    /**
	        * \brief Mark the bounding box as invalid.
	        *
	        * This operation sets the components of the minimum
	        * and maximum position to \f$\infty\f$ and \f$-\infty\f$,
	        * respectively.
	        */
	    inline void reset() {
		    min = PointType( std::numeric_limits<Scalar>::infinity());
		    max = PointType(-std::numeric_limits<Scalar>::infinity());
	    }

	    /// Calculate the n-dimensional volume of the bounding box
	    inline Scalar getVolume() const {
		    VectorType diff = max-min;
		    Scalar result = diff[0];
		    for (int i=1; i<PointType::dim; ++i)
			    result *= diff[i];
		    return result;
	    }

	    /// Calculate the n-1 dimensional volume of the boundary
	    inline Float getSurfaceArea() const {
		    VectorType d = max - min;
		    Float result = 0.0f;
		    for (int i=0; i<PointType::dim; ++i) {
			    Float term = 1.0f;
			    for (int j=0; j<PointType::dim; ++j) {
				    if (i == j)
					    continue;
				    term *= d[j];
			    }
			    result += term;
		    }
		    return 2.0f * result;
	    }

	    /// Return the center point
	    inline PointType getCenter() const {
		    return (max + min) * (Scalar) 0.5;
	    }

	    /// Check whether a point lies on or inside the bounding box
	    inline bool contains(const PointType &vec) const {
		    for (int i=0; i<PointType::dim; ++i)
			    if (vec[i] < min[i] || vec[i] > max[i])
				    return false;
		    return true;
	    }

	    /// Check whether a given bounding box is contained within this one
	    inline bool contains(const TAABB &aabb) const {
		    if (!isValid())
			    return false;
		    for (int i=0; i<PointType::dim; ++i)
			    if (aabb.min[i] < min[i] || aabb.max[i] > max[i])
				    return false;
		    return true;
	    }

	    /// Axis-aligned bounding box overlap test
	    inline bool overlaps(const TAABB &aabb) const {
		    for (int i=0; i<PointType::dim; ++i)
			    if (max[i] < aabb.min[i] || min[i] > aabb.max[i])
				    return false;
		    return true;
	    }

	    /// Expand the bounding box to contain another point
	    inline void expandBy(const PointType &p) {
		    for (int i=0; i<PointType::dim; ++i) {
			    min[i] = std::min(min[i], p[i]);
			    max[i] = std::max(max[i], p[i]);
		    }
	    }

	    /// Expand the bounding box to contain another bounding box
	    inline void expandBy(const TAABB &aabb) {
		    for (int i=0; i<PointType::dim; ++i) {
			    min[i] = std::min(min[i], aabb.min[i]);
			    max[i] = std::max(max[i], aabb.max[i]);
		    }
	    }

	    /// Calculate the squared point-AABB distance
	    inline Scalar squaredDistanceTo(const PointType &p) const {
		    Scalar result = 0;
		    for (int i=0; i<PointType::dim; ++i) {
			    Scalar value = 0;
			    if (p[i] < min[i])
				    value = min[i] - p[i];
			    else if (p[i] > max[i])
				    value = p[i] - max[i];
			    result += value*value;
		    }
		    return result;
	    }

	    /// Calculate the point-AABB distance
	    inline Scalar distanceTo(const PointType &p) const {
		    return std::sqrt(squaredDistanceTo(p));
	    }

	    /// Calculate the minimum squared AABB-AABB distance
	    inline Scalar squaredDistanceTo(const TAABB &aabb) const {
		    Scalar result = 0;

		    for (int i=0; i<PointType::dim; ++i) {
			    Scalar value = 0;
			    if (aabb.max[i] < min[i])
				    value = min[i] - aabb.max[i];
			    else if (aabb.min[i] > max[i])
				    value = aabb.min[i] - max[i];
			    result += value*value;
		    }
		    return result;
	    }

	    /// Calculate the minimum AABB-AABB distance
	    inline Scalar distanceTo(const TAABB &aabb) const {
		    return std::sqrt(squaredDistanceTo(aabb));
	    }

	    /// Return whether this bounding box is valid
	    inline bool isValid() const {
		    for (int i=0; i<PointType::dim; ++i)
			    if (max[i] < min[i])
				    return false;
		    return true;
	    }

	    /**
	        * \brief Return whether or not this bounding box
	        * covers anything at all.
	        *
	        * A bounding box which only covers a single point
	        * is considered nonempty.
	        */
	    inline bool isEmpty() const {
		    for (int i=0; i<PointType::dim; ++i) {
			    if (max[i] > min[i])
				    return false;
		    }
		    return true;
	    }

	    /// Return the axis index with the largest associated side length
	    inline int getLargestAxis() const {
		    VectorType d = max - min;
		    int largest = 0;

		    for (int i=1; i<PointType::dim; ++i)
			    if (d[i] > d[largest])
				    largest = i;
		    return largest;
	    }

	    /// Return the axis index with the shortest associated side length
	    inline int getShortestAxis() const {
		    VectorType d = max - min;
		    int shortest = 0;

		    for (int i=1; i<PointType::dim; ++i)
			    if (d[i] < d[shortest])
				    shortest = i;
		    return shortest;
	    }

	    /**
	        * \brief Calculate the bounding box extents
	        * \return max-min
	        */
	    inline VectorType getExtents() const {
		    return max - min;
	    }

	    /// Return a string representation of the bounding box
	    std::string toString() const {
		    std::ostringstream oss;
		    oss << "AABB[";
		    if (!isValid()) {
			    oss << "invalid";
		    } else {
			    oss << "min=" << min.toString()
				    << ", max=" << max.toString();
		    }
		    oss	<< "]";
		    return oss.str();
	    }

	    PointType min; ///< Component-wise minimum
	    PointType max; ///< Component-wise maximum
    };

    struct DatasetRecord {
        DatasetRecord( const Importance::Particle & particle ) {
            dir     = particle.incidentDir;
            value   = particle.weight;
            n       = 1;
            invDist = 1.f / particle.distance;
            fa      = 0.f;
        }

        DatasetRecord() {
            value   = 0.f;
            dir     = Vector3( 0.f );
            n       = 1;

            tmpHarmDist = 0.f;
            invDist = 1.f;
            fa      = 0.f;
        }

        void reset() {
            value   = 0.f;
            dir     = Vector3( 0.f );
            n       = 0;

            tmpHarmDist = 0.f;
            invDist = 1.f;
            fa      = 0.f;
        }

        void add( const Importance::Particle & particle ) {
            //dir    = ( (Float) n * dir + particle.direction ) / (Float) ( n + 1 );
            ++n;        
            //value  = (Float) n;
            
            dir     = ( value * dir + particle.weight * particle.incidentDir ) / ( value + particle.weight );
            invDist = value / invDist + particle.weight * particle.distance;            
            value   += particle.weight;
            invDist = value / invDist;                        
        }

        void merge( const DatasetRecord & rec ) {
            //dir    = ( (Float) n * dir + (Float) rec.n * rec.dir ) / (Float) ( n + rec.n );
            n += rec.n;                    
            //value  = (Float) n;

            dir     = ( value * dir + rec.value * rec.dir ) / ( value + rec.value );
            invDist = value / invDist + rec.value / rec.invDist;
            value   += rec.value;                        
            invDist = value / invDist;
        }


        //////////////////////////////////////////////////////////////////////////
        /// Adaption for Stepwise E-M
        //////////////////////////////////////////////////////////////////////////


        IMPORTANCE_INLINE Float getValue() const {
            return value;
        }

        IMPORTANCE_INLINE const Vector3 & getPoint() const {
            return dir;
        }


        IMPORTANCE_INLINE Float getInverseDistance() const {
            return invDist;
        }

        IMPORTANCE_INLINE Float getDistance() const { return 1.f / invDist; }

        //////////////////////////////////////////////////////////////////////////
        /// Data
        //////////////////////////////////////////////////////////////////////////

        /* */
        Float value;

        /* Direction or mean direction (in case of Verbeek) */
        Vector3 dir;

        /* Number of photons in the cell (for Verbeek otherwise its 1)  */
        unsigned int n;

        /* Temporary storage for harmonic distance computation */
        Float tmpHarmDist;

        Float invDist;

        /* Needed for a greedy version of Verbeek algorithm */
        Float fa;

        TAABB<Vector2> aabb;
    };
}