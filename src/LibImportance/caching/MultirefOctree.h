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
#include "MultirefOctreeAux.h"
#include <map>
#include "../shared/BoundingBox.h"

namespace Importance {

/// \brief Acceleration structure for searching of records influencing a point. Based 
///        on multiple-reference octree described by Krivanek.
class MultirefOctree {
    typedef BallRecord  Record;
    typedef MultirefOctreeNode<int> Node;
protected:

    template<class TOutIterator>
    void createChildBboxes(const BoundingBox3& input, TOutIterator& output) const {
        const Vector3 center = input.getCenter();
        for(int i = 0; i < 8; ++i) {
            output[i] = input;
            if ((i & 4) != 0) {
                output[i].point1.x = center.x;
            } else {
                output[i].point2.x = center.x;
            }
            if ((i & 2) != 0) {
                output[i].point1.y = center.y;
            } else {
                output[i].point2.y = center.y;
            }
            if ((i & 1) != 0) {
                output[i].point1.z = center.z;
            } else {
                output[i].point2.z = center.z;
            }
        }
    }

    /// \brief All nodes of the tree
    IStack<Node> nodes;

    /// \brief All points stored in the tree
    IStack<Record> records;

    /// \brief BoundingBox of the entire tree
    BoundingBox3 entireBbox;

    /// \brief maximal number of records in a leaf
    static const int leafSize = 15;

    /// \brief maximal depth of the tree
    int maxDepth;

    /// \brief Tests if it is feasible to insert given sphere into leaves of Octree node with 
    ///        given BoundingBox or if it is better to insert it into the inner node
    /// \param sphere Vector3 with its validity radius to insert
    /// \param bbox bounding box of the parent node
    /// \return true if the point should be inserted into current node and not into its 
    ///         children
    bool nonSplittable( const BoundingSphere& sphere, const BoundingBox3& bbox ) const {
        Vector3 bboxCenter = bbox.getCenter();
        Vector3 tested((bboxCenter.x > sphere.center.x) ? bbox.point2.x : bbox.point1.x,
                       (bboxCenter.y > sphere.center.y) ? bbox.point2.y : bbox.point1.y,
                       (bboxCenter.z > sphere.center.z) ? bbox.point2.z : bbox.point1.z);

        return (tested - sphere.center).square() <= 2.f*sphere.getRadiusSquared();
    }

    /// \brief Recursive procedure for building the tree. One call splits one node and 
    ///        recursively creates its children, or sets the node as leaf.
    /// \param nodeIndex index of the node to be split. The node is already allocated and 
    ///        contains all its points before the call.
    /// \param bBox Bounding Box of the current node
    /// \param depth depth of the node to be created
    void iterate( const int nodeIndex, const BoundingBox3& bBox, const int depth );

    /// \brief Adds single record to the tree, recursively propagating it into multiple 
    ///        subtrees if necessary
    /// \param nodeId index of tree node to add to (in the nodes stack)
    /// \param recordId index of the record to add (in the records stack)
    /// \param bbox bounding box of the node nodeId
    /// \param depth depth of the node nodeId
    void addRecordRecursively( const int nodeId, const int recordId, const BoundingBox3& bbox, 
                               const int depth );

    /// \brief updates maximal depth of this tree depending on number of records stored
    IMPORTANCE_INLINE void updateMaxDepth() {
        this->maxDepth = int(3+1.3f*Ff::log(float(records.size()/float(leafSize)), 8.f));
    }

public:

    /// \brief IIterator over query result set, which implements the whole search
    class IIterator {
    protected:

        /// \brief Tree that was queried
        const MultirefOctree* tree;

        /// \brief query point
        Vector3 query;

        /// \brief pointer to actual node during the iteration
        const Node* actual;

        /// \brief position in actual node records
        int position;

        /// \brief true if there are still more nodes to report
        bool isntFinished;

        /// \brief moves the iterator to next point in the query result
        IMPORTANCE_INLINE void next() {
            ++position;
            if(position >= actual->records.size()) {
                position = 0;
                do {
                    if(actual->isLeaf()) {
                        this->isntFinished = false;
                        return;
                    }
                    Vector3 actualCenter = actual->center;
                    int offset = 0;
                    offset += ((actualCenter.x < query.x) ? 4 : 0);
                    offset += ((actualCenter.y < query.y) ? 2 : 0);
                    offset += ((actualCenter.z < query.z) ? 1 : 0);
                    actual = & this->tree->nodes[actual->getChildPtr() + offset];        
                } while (actual->records.size() == 0);
            }
        }

    public:

        /// \brief Constructs new instance and starts the search
        /// \param tree where to search
        /// \param query location of the query
        /// \param kNN how many nearest neighbours to report (when using kNN search)
        IMPORTANCE_INLINE IIterator(const MultirefOctree* tree, const Vector3 query) {
            this->tree = tree;
            this->query = query;
            this->actual = &tree->nodes[0];
            this->position = -1;
            this->isntFinished = true;
            next();
        }

        /// \brief Returns true if there are any more results for the query.
        /// \return True if there are any more results for the query
        IMPORTANCE_INLINE bool hasNext() const {
            return this->isntFinished;
        }

        /// \brief returns next point from the result set and moves the iterator to next one
        /// \return index of current point from the result set
        IMPORTANCE_INLINE int getNext() {
            IMPORTANCE_ASSERT(hasNext());
            const int retVal = actual->records[position];
            next();
            return retVal;
        }
    };

    /// \brief Deallocates the tree
    void clear() {
        this->nodes.dealloc();
        this->records.dealloc();
    }

    /// \brief builds the tree over the specified set of records
    /// \param records records which will be stored in the tree
    /// \param bBox desired bounding box of the entire tree
    /// \param radiusQuocient constant by which to multiply each record R_i radius to obtain 
    ///                       their the validity radius
    template<class TRecord>
    void build(const IStack<TRecord>& records, const BoundingBox3& bBox, const float validityRadiusMult) {
        this->clear();
        this->entireBbox = bBox.getEpsilonEnlarged(1e-3f);
        this->nodes.reserve(1000);
        this->nodes.push(Node());
        this->nodes[0].records.reserve(records.size());

        for(int i = 0; i < records.size(); ++i) {
            this->records.push(Record(records[i]->position, records[i]->validRadius()*validityRadiusMult));
            this->nodes[0].records.push(i);
        }
        this->updateMaxDepth();
        iterate(0, entireBbox, 0);
    }

    /// \brief adds single record to the tree
    /// \param position position of the record
    /// \param normal normal of the record
    /// \param radius influence radius of the record
    void addRecord(const Vector3 position, const float radius, const float validityRadiusMult);

};
}