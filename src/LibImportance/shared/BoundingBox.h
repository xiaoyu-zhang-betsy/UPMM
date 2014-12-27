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


#ifndef ___BOUNDING_BOX___
#define ___BOUNDING_BOX___

#include "Vector3.h"

namespace Importance {

    /// \brief Axis aligned bounding box (AABB)
    class BoundingBox3 {
    public:

        /// \brief point with smallest coordinates in every direction
        Vector3 point1;        

        /// \brief point with biggest coordinates in every direction
        Vector3 point2;

        /// \brief Constructs default empty AABB
        IMPORTANCE_INLINE BoundingBox3() : point1(Vector3(INFINITY)), point2(Vector3(-INFINITY)) {}

        /// \brief Constructs AABB tightly enclosing 2 input AABBs. They do not need to be 
        ///        ordered in any way
        IMPORTANCE_INLINE BoundingBox3(const BoundingBox3& first, const BoundingBox3& second ) :
            point1(Vector3::min(first.point1, second.point1)), point2(Vector3::max(first.point2, second.point2)) {}


        /// \brief Constructs AABB tightly enclosing 2 arbitrary points. They do not need to be
        ///        ordered in any way
        IMPORTANCE_INLINE BoundingBox3( const Vector3 point1, const Vector3 point2 ) : 
            point1(Vector3::min(point1, point2)), point2(Vector3::max(point1, point2)) {}


        /// \brief Adds another BoundingBox to the current one, creating one that overlaps 
        ///        them both
        /// \param other BoundingBox to add
        /// \return self
        IMPORTANCE_INLINE BoundingBox3& operator+=( const BoundingBox3& other ) {
            this->point1 = Vector3::min(point1, other.point1);
            this->point2 = Vector3::max(point2, other.point2);
            return *this;
        }

        /// \brief adds another Vector2 to this box, enlarging this box as necessary 
        /// \param point point to add
        /// \return self
        IMPORTANCE_INLINE BoundingBox3& operator+=( const Vector3 point ) {
            this->point1 = Vector3::min(this->point1, point);
            this->point2 = Vector3::max(this->point2, point);
            return *this;
        }

        /// \brief Returns bounding box tightly enclosing this and second bounding box
        /// \param other bounding box to add
        /// \return bounding box tightly enclosing this and second bounding box
        IMPORTANCE_INLINE BoundingBox3 operator+( const BoundingBox3& other ) const {
            return BoundingBox3(Vector3::min(point1, other.point1), Vector3::max(point2, other.point2));
        }

        /// \brief Tests whether a BoundingBox is entirely contained in this one.
        /// \param other BoundingBox to test for containment
        /// \return true if other BoundingBox is entirely contained in this one
        IMPORTANCE_INLINE bool contains( const BoundingBox3& other ) const {        //TODO sse verze
            for(int i = 0; i < 3; ++i) {
                if(point1[i] > other.point1[i] || point2[i] < other.point2[i]) {
                    return false;
                }
            }
            return true;
        }
        /// \brief tells if a given point is inside BoundingBox
        /// \param point point to evaluate
        /// \return true if point is inside BoundingBox
        IMPORTANCE_INLINE bool contains(const Vector3 point) const {        //TODO SSE verze
            for(int i = 0; i < 3; ++i) {
                if(point[i] > point2[i] || point[i] < point1[i]) {
                    return false;
                }
            }
            return true;
        }

        /// \brief Returns intersection of this bounding box with another. Result should be 
        ///        tested if it is empty before using
        /// \param other AABB to intersect
        /// \return Intersection of this and other bounding box
        IMPORTANCE_INLINE BoundingBox3 getIntersection( const BoundingBox3& other ) const {
            BoundingBox3 result;
            result.point1 = Vector3::max(this->point1, other.point1);
            result.point2 = Vector3::min(this->point2, other.point2);
            return result;
        }



        /// \brief returns size of this BoundingBox
        /// \return vector with size of this BoundingBox in every direction
        IMPORTANCE_INLINE Vector3 size() const {
            return Vector3(point2 - point1);
        }

        /// \brief returns center point (centroid) of this BoundingBox
        /// \return center of this BoundingBox
        IMPORTANCE_INLINE Vector3 getCenter() const {
            return Vector3( (point2+point1) * 0.5f );
        }

    
        /// \brief compares this bounding box with another. Differences smaller than EPSILON 
        ///        are tolerated
        /// \param other BoundingBox to compare with
        /// \return true if vertices of the two BoundingBoxes differ by at most epsilon
        IMPORTANCE_INLINE bool operator==(const BoundingBox3& other) const {
            return this->point1 == other.point1 && this->point2 == other.point2;
        }

        /// \brief Sets this BoundingBox to be empty. (from infinity to -infinity). Such 
        ///        BoundingBox is only a dummy and has undefined certain operations, such as 
        ///        getVoume, getArea, etc.
        IMPORTANCE_INLINE void setEmpty() {
            this->point1 = Vector3(INFINITY);
            this->point2 = Vector3(-INFINITY);        
        }

        /// \brief Tests if this BoundingBox is empty
        /// \return true if this BoundingBox is empty
        IMPORTANCE_INLINE bool isEmpty() const {   //TODO sse verze
            for(int i = 0; i < 3; ++i) {
                if(point1[i] > point2[i]) {
                    return true;
                }
            }
            return false;
        }

        /// \brief returns squared shortest distance from this BoundingBox to a point
        /// \param vect point to which to measure the distance
        /// \return squared shortest distance from this BoundingBox to a point
        IMPORTANCE_INLINE Float distanceSquared(const Vector3 vect) const {
            return Vector3::max(Vector3::max(point1 - vect, vect - point2), Vector3(0.f)).square();
        }

        /// \brief returns squared largest distance from this BoundingBox to a point
        /// \param vect point to which to measure the distance
        /// \return squared largest distance from this BoundingBox to a point
        IMPORTANCE_INLINE Float distanceSquaredMax(const Vector3 vect) const {
            return Vector3::max(vect - point1, point2 - vect).square();
        }

        /// \brief Returns copy of this bounding box enlarged in each dimension by AABB_EPSILON
        /// \return copy of this bounding box enlarged in each dimension by AABB_EPSILON
        IMPORTANCE_INLINE BoundingBox3 getEpsilonEnlarged(const Float factor) const {
            Vector3 center = this->getCenter();
            Vector3 size = this->size();
            size *= (1.f+factor)*0.5f;
            size += Vector3(factor);
            return BoundingBox3(center-size, center+size);
        }

        /// \brief Returns smallest axis-aligned cube enclosing this bounding box
        /// \return smallest axis-aligned cube enclosing this bounding box
        IMPORTANCE_INLINE BoundingBox3 getCube() const {
            const Vector3 center = this->getCenter();
            const Float size = this->size().max()*0.5f;
            return BoundingBox3(center-Vector3(size), center+Vector3(size));
        }

        IMPORTANCE_INLINE bool isReal() const {
            return point1.isReal() && point2.isReal();
        }
    };
}

#endif