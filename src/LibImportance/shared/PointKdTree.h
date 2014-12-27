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


#ifndef ___POINT_KD_TREE__
#define ___POINT_KD_TREE__

#pragma warning(push)
#pragma warning(disable:4127)

#include "PointKdTreeNode.h"
#include "PointKdTreeShared.h"
#include "BoundingBox.h"
#include "StaticHeap.h"

namespace Importance {

    struct Kd3PositionTraits {
        typedef PointKdTreeNode Node;
        typedef Vector3         Point;
        template<class TData>
        static IMPORTANCE_INLINE Point extractPosition(const TData& data) {
            return data.position;
        }
        static IMPORTANCE_INLINE int dimension() {
            return 3;
        }
    };

    struct Kd3VectorTraits {
        typedef PointKdTreeNode Node;
        typedef Vector3         Point;
        template<class TData>
        static IMPORTANCE_INLINE Point extractPosition(const TData& data) {
            return data;
        }
        static IMPORTANCE_INLINE int dimension() {
            return 3;
        }
    };

    struct Kd2DirTraits {
        typedef PointKdTreeNode Node;
        typedef Vector2         Point;
        template<class TData>
        static IMPORTANCE_INLINE Point extractPosition(const TData& data) {
            return data.dir;
        }
        static IMPORTANCE_INLINE int dimension() {
            return 3;
        }
    };


    /// \brief Static kD-tree data structure for quick querying of point data. Stores 3D point
    ///        data and supports range-limited kNN queries. Uses longest side cut with sliding 
    ///        midpoint, tracking nodes and priority node searching. 
    ///        
    /// Its time complexity for kNN query is O(log n + k*log k). It is built in O(n log n) time.
    /// Because it uses packed nodes (PointKdTreeNode), its maximum node size is 127 and can have
    /// only 2^25 leaves
    template<typename Traits>
    class PointKdTree {
        typedef typename Traits::Node  Node;
        typedef typename Traits::Point Point;
        typedef TreeRecordBase<Point>  Record;
        typedef OpenKdNodeBase<Point>  OpenNode;
        typedef BoundingBox3     BBox;

    protected:

        /// \brief Bounding box of the entire scene
        BBox entireBBox;


        /// \brief Storage for all points in the tree
        IStack<Record> data;

        /// \brief Storage of all the nodes of the tree
        IStack<Node> nodes;

        /// \brief Bounding box stored for each leaf of the tree, used to speed up search by 
        ///        possibly rejecting leaf nodes which bounding box is further from the query 
        ///        point than k-th nearest neighbor already found
        IStack<BBox> leafBBoxes;


        /// \brief Maximum number of points stored in each leaf
        int leafSize;


        /// \brief Creates a leaf node from given range of points and links it to the parent node
        /// \param parent index of the parent node
        /// \param from first index of point from the data array to be included in this leaf
        /// \param to index past the last point from the data array to be included
        void createLeafNode(const unsigned int parent, const unsigned int from, const unsigned int to) {
            IMPORTANCE_ASSERT((to-from > 0 || this->data.size() == 0) && to-from <= Node::MAXIMUM_LEAF_SIZE);

            Node leaf;
            leaf.setIsLeaf();
            leaf.setLeafStart(from);
            leaf.setLeafSize(to-from);
            leaf.setLeafNumber(unsigned int(leafBBoxes.size()));

            BBox tightBbox;
            for(unsigned int i = from; i < to; ++i) {
                tightBbox += this->data[i].position;
            }
            leafBBoxes.push(tightBbox);

            unsigned int index = unsigned int(nodes.size());
            nodes.push(leaf);

            if(index - parent == 1 || parent == -1) {
                // left child is implicit
                return;
            } else {
                IMPORTANCE_ASSERT(!nodes[parent].hasRightChild());
                nodes[parent].setRightChild(index);
            }
        }

        /// \brief Recursively creates a subtree from given range of points and links it to the 
        ///        parent node
        /// \param parent parent node of this subtree
        /// \param from first index of point from the data array to be included in this subtree
        /// \param to index past the last point from the data array to be included
        /// \param box Bounding box associated with this subtree
        void iterate(const unsigned int parent, const int from, const int to, const BBox& box) {
            for(int i = from; i < to; ++i) {
                IMPORTANCE_ASSERT(box.contains(data[i].position));
            }
            int splitDim = -1;
            bool terminate = false;
            Point cellSize = box.size();

            if(to - from > leafSize) {
                //find largest axis in which are the data non-singular
                for(int i = 0; i < Traits::dimension(); ++i ) {
                    splitDim = cellSize.argMax();
                    bool dimensionNonSingular = false;
                    for(int j = from; j < to; ++j) {
                        if(data[to-1].position[splitDim] != data[j].position[splitDim]) {
                            dimensionNonSingular = true;
                            break;
                        }
                    }
                    if(dimensionNonSingular) {
                        break;
                    } else {
                        cellSize[splitDim] = 0.f;
                        if(i == Traits::dimension()-1) {
                            terminate = true;
                        }
                    }

                }
            } else {
                terminate = true;
            }

            if(terminate && (to-from) < Node::MAXIMUM_LEAF_SIZE) {
                createLeafNode(parent, from, to);
                return;
            }

            Float split = box.getCenter()[splitDim];
            int leftPtr = from;
            int rightPtr = to-1;

            BBox bboxL = box;
            BBox bboxR = box;


            // partition the data
            if(!terminate) {
                while(true) {
                    while(leftPtr < to && data[leftPtr].position[splitDim] <= split) {
                        ++leftPtr;
                    }
                    while(rightPtr >= from && data[rightPtr].position[splitDim] >= split ) {
                        --rightPtr;
                    }

                    if(leftPtr < rightPtr) {
                        std::swap(data[leftPtr], data[rightPtr]);
                    } else {
                        break;
                    }
                }

                // leftPtr = now first on right
                if( leftPtr == from ) {   // left is empty -> sliding midpoint
                    unsigned int minIndex = from;
                    split = data[from].position[splitDim];
                    for(int i = from+1; i < to; ++i) {
                        if(data[i].position[splitDim] < split) {
                            minIndex = i;
                            split = data[i].position[splitDim];
                        }
                    }
                    std::swap(data[from], data[minIndex]);

                    ++leftPtr;
                } else if( leftPtr == to ) {   // right is empty -> sliding midpoint
                    unsigned int maxIndex = from;
                    split = data[from].position[splitDim];
                    for(int i = from+1; i < to; ++i) {
                        if(data[i].position[splitDim] > split) {
                            maxIndex = i;
                            split = data[i].position[splitDim];
                        }
                    }
                    std::swap(data[to-1], data[maxIndex]);
                    --leftPtr;
                }

                bboxR.point1[splitDim] = bboxL.point2[splitDim] = split;

            } else {
                leftPtr = (to+from)/2;
                bboxR.point1[splitDim] = bboxL.point2[splitDim] = data[leftPtr].position[splitDim];
            }

            for(int i = from; i < leftPtr; ++i) {
                IMPORTANCE_ASSERT(bboxL.contains(data[i].position));
            }
            for(int i = leftPtr; i < to; ++i) {
                IMPORTANCE_ASSERT(bboxR.contains(data[i].position));
            }

            Node node;
            node.setSplitAxis(splitDim);
            node.setSplitPosition(float(split));
            unsigned int index = int(nodes.size());
            nodes.push(node);
            if(index - parent != 1) { // left is implicit
                IMPORTANCE_ASSERT(!nodes[parent].hasRightChild());
                nodes[parent].setRightChild(index);
            }

            iterate (index, from,    leftPtr, bboxL);
            iterate (index, leftPtr, to,      bboxR);
        }


    public:        
        class AlwaysTruePredicate {
        public:        
            bool operator()( const unsigned int & ) const { return true; }
        };

        /// \brief Deallocates the tree and frees the memory
        void dealloc() {
            this->data.dealloc();
            this->nodes.dealloc();
            this->leafBBoxes.dealloc();
            this->entireBBox.setEmpty();
        }

        /// \brief Builds the tree over given input elements. The tree will store the data in form 
        ///        of indices based on position of each element in the input array. The method is 
        ///        templated to work with Vector3 array or array of any data types which have
        ///        public Vector3 member position.
        /// \param items pointer to the first element of the continuous storage where the input 
        ///              data are located
        /// \param size number of elements in the input array
        /// \param leafSize Maximum number of elements to be stored in single leaf
        template<class TDataIterator>
        void build(const TDataIterator begin, const TDataIterator end, const int leafSize) {
            this->leafSize = leafSize;
            this->nodes.clear();
            this->leafBBoxes.clear();
            this->data.resize(int(end-begin));
            entireBBox.setEmpty();

        
            int ptr = 0;
            for(auto it = begin; it != end; it++) {
                this->data[ptr].data = ptr;
                this->data[ptr].position = Traits::extractPosition(*it);
                entireBBox += this->data[ptr].position;
                ptr++;
            }
            iterate(unsigned(-1), 0, int(data.size()), entireBBox);
        }

        /// \brief Performs range limited kNN query in the tree and reports indices of found 
        ///        points as well as squared distances from the query to each result. Photon returned at the zero 
        ///        index in the array is the furthest one of all returned
        /// \param origin Query point position
        /// \param range maximum range of the query (query is therefore limited to a ball volume)
        /// \param _results where to write results
        /// \param maxCount how many nearest neighbors to report
        /// \param _found Optional parameter telling how many nearest neighbors are already stored 
        ///               in the results array. This parameter allows to make incremental kNN 
        ///               queries with increasing k
        /// \return number of nearest neighbors reported
        template<class TPredicate>
        int rangeQuery(const Point origin, const Float range, KdQueryResult* _results, const int maxCount, const TPredicate & allowParticle, int _found=0) const {
            StaticHeap<OpenNode, 1500> openNodes;
            ImExternalHeap<KdQueryResult, KdQueryResultComparator> heap(_results, maxCount, _found);
            Float maxDistSqr = range*range;  
            Point distSqr = Point::max(Point(0.f), Point::max(entireBBox.point1-origin, origin-entireBBox.point2));                 
            distSqr *= distSqr;

            OpenNode actual(0, distSqr.l1Norm(), distSqr);
            Node actualNode = nodes[actual.node];

            while(actual.distSqr < maxDistSqr) {
                if(actualNode.isLeaf()) {
                    if(leafBBoxes[actualNode.getLeafNumber()].distanceSquared(origin) < maxDistSqr) {

                        const unsigned int to = actualNode.getLeafStart() + actualNode.getLeafSize();
                        for(unsigned int i = actualNode.getLeafStart(); i < to; ++i) {
                            // insert all points in the leaf to max-heap in results array
                            const Float distSquared = (origin-this->data[i].position).square();
                            if(distSquared < maxDistSqr && allowParticle( this->data[i].data )) {
                                const KdQueryResult enqueued(data[i].data, distSquared);
                                if(heap.isFull()) { // the heap is full, overwrite first and repairTop
                                    heap.replaceTop(enqueued);
                                    maxDistSqr = std::min(heap.peek().distSqr, maxDistSqr);
                                } else {                // heap insert
                                    heap.insert(enqueued);
                                    if(heap.isFull()) {
                                        maxDistSqr = std::min(heap.peek().distSqr, maxDistSqr);
                                    }
                                }
                            }
                        }
                    }
                    if(openNodes.isEmpty()) {
                        break;
                    }
                    actual = openNodes.removeTop();
                    actualNode = nodes[actual.node];

                } else {    // inner node
                    const int splitAxis = actualNode.getSplitAxis();
                    const Float splitPosition = actualNode.getSplitPosition();

                    if(origin[splitAxis] < splitPosition) {  //left is closer
                        OpenNode right = actual;
                        actual = OpenNode(actual.node+1, actual.distSqr, actual.distancesSqr);
                        actualNode = nodes[actual.node];
                        const Float newOffset = splitPosition - origin[splitAxis];
                        right.distSqr -= right.distancesSqr[splitAxis];
                        right.distancesSqr[splitAxis] = newOffset*newOffset;
                        right.distSqr += right.distancesSqr[splitAxis];
                        if(right.distSqr < maxDistSqr) {
                            right.node = nodes[right.node].getRightChild();
                            openNodes.insert(right);
                        }
                    } else {                                 //right is closer
                        OpenNode left = actual;
                        actual = OpenNode(actualNode.getRightChild(), actual.distSqr, actual.distancesSqr);
                        actualNode = nodes[actual.node];
                        const Float newOffset = splitPosition - origin[splitAxis];
                        left.distSqr -= left.distancesSqr[splitAxis];
                        left.distancesSqr[splitAxis] = newOffset*newOffset;
                        left.distSqr += left.distancesSqr[splitAxis];
                        if(left.distSqr < maxDistSqr) {
                            left.node = left.node+1;
                            openNodes.insert(left);
                        }
                    }
                }
            }
            return heap.size();
        }
        
        IMPORTANCE_INLINE int rangeQuery(const Point origin, const Float range, KdQueryResult* _results, const int maxCount, int _found=0) const {
            return rangeQuery(origin, range, _results, maxCount, AlwaysTruePredicate(), _found );
        }
    };
}
#pragma warning(pop)

#endif
