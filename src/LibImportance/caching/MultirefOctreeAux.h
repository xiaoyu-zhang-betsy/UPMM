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
#include "../shared/Config.h"
#include "../shared/Vector3.h"
#include <vector>
#include "../shared/BoundingBox.h"
#include "../shared/Stack.h"

namespace Importance {

/// \brief Single node in an multiple-reference octree. Because there is no distinction 
///        between the node can both have children and contain records
template<class TRecord>
class MultirefOctreeNode {

protected:

    /// \brief index of first of child nodes, or -1 if this is leaf
    int childPtr;

public:
    /// \brief Center point of this node, to speed up traversal
    Vector3 center;

    /// \brief indices of all records stored in this node
    IStack<TRecord> records;

    /// \brief sets whether this is child node or not
    /// \param leaf true if this node is to be set a leaf
    IMPORTANCE_INLINE void setLeaf(const bool leaf) {
        if(leaf == true) {
            childPtr = -1;  // -1 is leaf
        } else {
            // if this was leaf, we need to set childPtr to something else than -1
            if(childPtr == - 1) {
                childPtr = -2;
            }
        }
    }

    /// \brief returns true if this is leaf node
    /// \return true if this is leaf node
    IMPORTANCE_INLINE bool isLeaf() const {
        return childPtr == -1;
    }

    /// \brief sets index of the first child node. Node must be inner
    /// \param ptr child node index to set
    IMPORTANCE_INLINE void setChildPtr(const int ptr) {
        IMPORTANCE_ASSERT(!isLeaf());
        this->childPtr = ptr;
    }

    /// \brief Returns index of the first child node. Node must be inner
    /// \return index of the first child node
    IMPORTANCE_INLINE int getChildPtr() const {
        IMPORTANCE_ASSERT(!isLeaf());
        return this->childPtr;
    }
};



/// \brief Single record in the multiple-reference IC Octree. Used because records in the 
///        octree can have different validity radius than in the IC
struct BallRecord {
protected:
    /// \brief record position in space
    Vector3 _position;
    
public:
    IMPORTANCE_INLINE BallRecord() { }

    IMPORTANCE_INLINE BallRecord(const Vector3 position, const float validity) {
        this->_position = position;
        _position.extraData = validity;
    }
    IMPORTANCE_INLINE Vector3 position() const {
        return _position;
    }
    IMPORTANCE_INLINE float radius() const {
        return _position.extraData;
    }
};

/// \brief Represents 3D bounding sphere.
class BoundingSphere {
protected:

    /// \brief radius of this sphere
    float radius;

    /// \brief radius of this sphere squared
    float radiusSquared;

public:

    /// \brief center point of sphere
    Vector3 center;    


    /// \brief Constructs default BoundingSphere
    IMPORTANCE_INLINE BoundingSphere() {
#ifdef LIBIMP_DEBUG
        radiusSquared = radius = NAN;
#endif
    }

    /// \brief Constructs BoundingSphere with given parameters
    /// \param center center point of sphere
    /// \param radius radius of sphere. Should be positive, else this class is not guaranteed 
    ///               to work
    IMPORTANCE_INLINE BoundingSphere( const Vector3 center, const float radius ) : radius(radius), radiusSquared(radius*radius), center(center) {
        IMPORTANCE_ASSERT(radius >= 0);        
    }

    /// \brief getter for radius
    /// \return radius of this sphere
    IMPORTANCE_INLINE float getRadius() const {
        return radius;
    }

    /// \brief getter for radius squared
    /// \return radius^2
    IMPORTANCE_INLINE float getRadiusSquared() const {
        return radiusSquared;
    }

    /// \brief Sets radius squared. Radius is computed automatically.
    /// \param radiusSquared radius^2 to set
    IMPORTANCE_INLINE void setRadiusSquared(const float radiusSquared) {
        this->radiusSquared = radiusSquared;
        this->radius = sqrt(radiusSquared);
    }

    /// \brief Find out whether a given Vector2 lies within (or on the border of) this sphere
    /// \param point point to inspect
    /// \param distanceSquared output parameter, where to write distance^2 from point to center
    /// \return true if point is within this sphere
    IMPORTANCE_INLINE bool intersect(const Vector3 point, float& distanceSquared) const {
        distanceSquared = (point - this->center).square();
        return distanceSquared <= radiusSquared;
    }

    /// \brief Determine whether this sphere intersects a given AABB
    /// \param box bounding box to inspect
    /// \return true if box intersects this sphere
    IMPORTANCE_INLINE bool intersect(const BoundingBox3& box) const {
        float dummy;
        return intersect(box, dummy);
    }

    //TODO pokud to nekde pouzivam, tak prepsat na SSE
    /// \brief Determine whether this sphere intersects a given AABB
    /// \param box bounding box to inspect
    /// \return true if box intersects this sphere
    IMPORTANCE_INLINE BoundingBox3 getIntersection(const BoundingBox3& box) const {
        BoundingBox3 sphereBox = this->getBbox();
        //BoundingBox3 TODObbox = sphereBox;
        for(int ax = 0; ax <= 2; ++ax) {
            float cutRadius;
            if(box.point1[ax] < this->center[ax]) {
                if(this->center[ax] < box.point2[ax]) { //je uprostred
                    cutRadius = getRadius();
                } else {  //je napravo
                    const float diff = this->center[ax] - box.point2[ax];
                    IMPORTANCE_ASSERT(diff >= 0);
                    if(diff > this->getRadius()) {
                        return BoundingBox3();
                    } else {
                        cutRadius = Ff::sqrtFast(getRadiusSquared() - diff*diff);
                    }
                }
            } else {    // je nalevo
                IMPORTANCE_ASSERT(this->center[ax] < box.point2[ax]);
                const float diff = box.point1[ax] - this->center[ax];
                IMPORTANCE_ASSERT(diff >= 0);
                if(diff > this->getRadius()) {
                    return BoundingBox3();
                } else {
                    cutRadius = Ff::sqrtFast(getRadiusSquared() - diff*diff);
                }
            }
            IMPORTANCE_ASSERT(isReal(cutRadius) && cutRadius >= 0 && cutRadius <= getRadius());
            const int axis1 = (ax+1)%3;
            const int axis2 = (ax+2)%3;
            sphereBox.point1[axis1] = std::max(sphereBox.point1[axis1], center[axis1]-cutRadius);
            sphereBox.point2[axis1] = std::min(sphereBox.point2[axis1], center[axis1]+cutRadius);
            sphereBox.point1[axis2] = std::max(sphereBox.point1[axis2], center[axis2]-cutRadius);
            sphereBox.point2[axis2] = std::min(sphereBox.point2[axis2], center[axis2]+cutRadius);
        }
        return sphereBox.getIntersection(box);
    }

    /// \brief Determine whether this sphere intersects a given AABB
    /// \param box bounding box to inspect
    /// \param distanceSquared output parameter, where to write minimal distance^2 from box 
    ///                        to center
    /// \return true if box intersects this sphere
    IMPORTANCE_INLINE bool intersect(const BoundingBox3& box, float& distanceSquared) const {  //TODO SSE
        distanceSquared = 0.f;
        for (int i = 0; i < 3; ++i) {
            float tempDist = std::max( std::max(box.point1[i]-center[i], center[i]-box.point2[i]), 0.f );
            distanceSquared += tempDist * tempDist;
        }
        return distanceSquared < radiusSquared;
    }

    /// \brief If a given Vector2 is within the sphere's radius, shrink the sphere so the point
    ///        lays on sphere surface
    /// \param p Vector2 which will now lay on the sphere surface
    IMPORTANCE_INLINE void shrink(const Vector3 p) {
        float distanceSquared = (this->center - p).square();
        if (distanceSquared < radiusSquared) {
            this->setRadiusSquared(radiusSquared);
        }
    }

    /// \brief Returns bounding box tightly enclosing this sphere
    /// \return bounding box tightly enclosing this sphere
    IMPORTANCE_INLINE BoundingBox3 getBbox() const {
        const Vector3 radiusVector(this->getRadius());
        return BoundingBox3(center-radiusVector, center+radiusVector);
    }

    IMPORTANCE_INLINE float getArea() const {
        return (4.f*PI) * radiusSquared;
    }


};

IMPORTANCE_INLINE BoundingBox3& operator+=( BoundingBox3& box, const BoundingSphere& sphere ) {
    return box += sphere.getBbox();
}


}