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


#ifndef ___DYNAMIC_KD_TREE__
#define ___DYNAMIC_KD_TREE__

#include "../shared/Vector3.h"
#include "../shared/PointKdTreeNode.h"
#include "../shared/PointKdTreeShared.h"
#include "../shared/BoundingBox.h"
#include <vector>

namespace Importance {

    struct DynamicKd3PositionTraits {
        typedef Dynamic3dTreeNode Node;
        typedef Vector3         Point;
        template<class TData>
        static IMPORTANCE_INLINE Point extractPosition(const TData& data) {
            return data.position;
        }
        static IMPORTANCE_INLINE int dimension() {
            return 3;
        }
    };

    struct DynamicKd3PtrPositionTraits {
        typedef Dynamic3dTreeNode Node;
        typedef Vector3         Point;
        template<class TData>
        static IMPORTANCE_INLINE Point extractPosition(const TData data) {
            return data->position;
        }
        static IMPORTANCE_INLINE int dimension() {
            return 3;
        }
    };

    struct DynamicKd3VectorTraits {
        typedef Dynamic3dTreeNode Node;
        typedef Vector3         Point;
        template<class TData>
        static IMPORTANCE_INLINE Point extractPosition(const TData& data) {
            return data;
        }
        static IMPORTANCE_INLINE int dimension() {
            return 3;
        }
    };


    template<typename Traits>
    class ImDynamicKdTree {
        typedef typename Traits::Node  Node;
        typedef typename Traits::Point Point;
        typedef TreeRecordBase<Point>  Record;
        typedef OpenKdNodeBase<Point>  OpenNode;
        typedef BoundingBox3           BBox;

        /// \brief Bounding box of the entire scene
        BBox entireBBox;


        /// \brief Structure of single dynamic kd tree leaf
        struct Leaf {

            /// \brief All points stored in this leaf
            IStack<Record> points;

            /// \brief Bounding box of this leaf
            BBox bbox;
        };

        /// \brief Storage of all leaves data of the tree
        IStack<Leaf> leaves;

        /// \brief Storage of all nodes of the tree (inner and leaves)
        IStack<Node> nodes;

        /// \brief Maximum number of points stored in each leaf
        int leafSize;

        /// \brief How many records are stored in the tree
        int stored;

        /// \brief Creates a leaf from given node and leaf data
        /// \param nodeIndex index of the node stored in nodes array to create the leaf from
        /// \param leafIndex index of leaf data in the leaves array to create the leaf from
        void createLeafNode( const int nodeIndex, const int leafIndex ) {
            Node& node = this->nodes[nodeIndex];
            node.setIsLeaf();
            node.setLeafData(leafIndex);
            BBox tightBbox;
            for(int i = 0; i < leaves[leafIndex].points.size(); ++i) {
                tightBbox += leaves[leafIndex].points[i].position;
            }
            leaves[leafIndex].bbox = tightBbox;
        }

        /// \brief Recursively creates a subtree from given leaf data and given node
        /// \param nodeIndex index of the node which to process
        /// \param leafIndex index of leaf data where the data for this subtree and their 
        ///                  bounding box is stored
        /// \param bbox bounding box of current node
        void iterate( const int nodeIndex, const int leafIndex, const BBox& bbox ) {
            IStack<Record>& leaf = this->leaves[leafIndex].points;

            IMPORTANCE_ASSERT(leaves.size() <= this->stored+2 && nodes.size() <= this->stored+2);

            int splitDim = -1;
            bool terminate = false;
            Point cellSize = bbox.size();

            if(leaf.size() > leafSize) {
                //find largest axis in which are the data non-singular
                for(int i = 0; i < Traits::dimension(); ++i ) {
                    splitDim = cellSize.argMax();
                    bool dimensionNonSingular = false;
                    for(int j = 0; j < leaf.size(); ++j) {
                        if(abs(leaf.back().position[splitDim]-leaf[j].position[splitDim]) > 1e-6f) {
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

            if(terminate ) {
                createLeafNode(nodeIndex, leafIndex);
                IMPORTANCE_ASSERT(nodes[nodeIndex].getLeafData() < this->leaves.size());
                return;
            }

            unsigned int rightIndex = unsigned int(leaves.size());
            Float split = bbox.getCenter()[splitDim];
            this->leaves.push(Leaf());
            IStack<Record>& leftLeaf = leaves[leafIndex].points;
            IStack<Record>& rightLeaf = leaves[rightIndex].points;
            IMPORTANCE_ASSERT(rightLeaf.size() == 0);

            int leftPtr = 0;
            // partition 
            for(int i = 0; i < leftLeaf.size(); ++i) {
                if(leftLeaf[i].position[splitDim] <= split) {
                    leftLeaf[leftPtr++] = leftLeaf[i];
                } else {
                    rightLeaf.push(leftLeaf[i]);
                }
            }
            leftLeaf.resize(leftPtr);

    #if 1
            // leftPtr = now first on right
            if( leftLeaf.size() == 0 ) {   // left is empty -> sliding midpoint
                int minIndex = -1;
                split = INFINITY;
                for(int i = 0; i < rightLeaf.size(); ++i) {
                    if(rightLeaf[i].position[splitDim] < split) {
                        minIndex = i;
                        split = rightLeaf[i].position[splitDim];
                    }
                }
                leftLeaf.push(rightLeaf[minIndex]);
                rightLeaf.remove(minIndex);
            } else if( rightLeaf.size() == 0 ) {   // right is empty -> sliding midpoint
                int maxIndex = -1;
                split = -INFINITY;
                for(int i = 0; i < leftLeaf.size(); ++i) {
                    if(leftLeaf[i].position[splitDim] > split) {
                        maxIndex = i;
                        split = leftLeaf[i].position[splitDim];
                    }
                }
                rightLeaf.push(leftLeaf[maxIndex]);
                leftLeaf.remove(maxIndex);
            }
    #endif

            BBox bboxL = bbox;
            BBox bboxR = bbox;
            bboxL.point2[splitDim] = split;
            bboxR.point1[splitDim] = split;
            IMPORTANCE_ASSERT(isReal(split));

            // push 2 children
            //nodes.resize(nodes.size()+2);
            nodes.push(Node());
            nodes.push(Node());
            Node& node = nodes[nodeIndex];
            node.setSplitAxis(splitDim);
            node.setSplitPosition(float(split));
            node.setLeftChild(int(nodes.size()-2));

            IMPORTANCE_ASSERT(nodes[nodeIndex].getLeftChild()+1 < nodes.size());
            const int nodeSize = int(nodes.size());
            iterate(nodeSize-2, leafIndex, bboxL);
            iterate(nodeSize-1, rightIndex, bboxR);
            IMPORTANCE_ASSERT(nodes[nodeIndex].getLeftChild()+1 < nodes.size());
        }

    public:
        
        /// \brief adds single sample to the tree
        /// \param sample the sample to add
        void add(const Vector3 position, const unsigned int /*data*/) {
            IMPORTANCE_ASSERT(this->entireBBox.contains(position));
            int nodeIndex = 0;
            BBox bbox = this->entireBBox;
            while (!nodes[nodeIndex].isLeaf()) {
                // find the leaf where the sample belongs
                const float splitPos = nodes[nodeIndex].getSplitPosition();
                const int splitAxis = nodes[nodeIndex].getSplitAxis();
                const bool left = position[splitAxis] < splitPos;
                if(left) {
                    bbox.point2[splitAxis] = splitPos;
                    nodeIndex = nodes[nodeIndex].getLeftChild();
                } else {
                    bbox.point1[splitAxis] = splitPos;
                    nodeIndex = nodes[nodeIndex].getLeftChild()+1;
                }
                IMPORTANCE_ASSERT(bbox.point1[splitAxis] <= bbox.point2[splitAxis]);
            }
            const Node& node = nodes[nodeIndex];
            Leaf& leaf = leaves[node.getLeafData()];

            leaf.points.push(Record(stored++, position));
            leaf.bbox += position;

            if(leaf.points.size() > this->leafSize) {
                iterate(nodeIndex, node.getLeafData(), bbox);
            }
        }

        /// \brief Clears contents of this tree and frees the memory
        void dealloc() {
            this->leaves.resize(0);
            this->nodes.resize(0);
        }

        /// \brief builds the tree over the specified set of samples
        /// \param samples samples which will be stored in the tree
        /// \param bBox desired bounding box of the entire tree
        template<class TItemIterator>
        void build(const TItemIterator begin, const TItemIterator end, const BoundingBox3& bBox, const int leafSize) {
            this->leafSize = leafSize;
            this->dealloc();
            this->leaves.push(Leaf());
            this->leaves[0].points.reserve(int(end-begin));
            int ptr = 0;
            for(auto it = begin; it != end; it++) {
                this->leaves[0].points.push(Record(ptr, Traits::extractPosition(*it)));
                ptr++;
            }


            this->entireBBox = bBox;

            entireBBox = entireBBox.getEpsilonEnlarged(0.001f);
            this->nodes.push(Node());
            this->stored = int(end-begin);
            iterate(0, 0, entireBBox);
        }

        int getClosest(const Vector3 position, const Float maxDist, const int kNN, KdQueryResult* output) const {
            IMPORTANCE_ASSERT(kNN > 0);
            if(nodes.size() == 0) {
                return 0;
            }
            StaticHeap<OpenNode, 5*100> openNodes;
            ImExternalHeap<KdQueryResult, KdQueryResultComparator> heap(output, kNN, 0);
            Float maxDistSqr = maxDist*maxDist;
            Point distSqr = Point::max(Point(0.f), Point::max(entireBBox.point1-position, position-entireBBox.point2));                 

            OpenNode actual = OpenNode(0, distSqr.l1Norm(), distSqr);
            Node actualNode = nodes[actual.node];

            const float EPS_FACTOR = 0.6f;

            while(actual.distSqr < EPS_FACTOR*maxDistSqr) {
                if(actualNode.isLeaf()) {
                    const Leaf& leaf = leaves[actualNode.getLeafData()];

                    for(int i = 0; i < leaf.points.size(); ++i) {
                        // insert all points in the leaf to max-heap in results array
                        const Float distSquared = (position - leaf.points[i].position).square();
                        if(distSquared < maxDistSqr) {
                            KdQueryResult inserted(leaf.points[i].data, distSquared);
                            if(heap.isFull()) { // the heap is full, overwrite first and repairTop
                                heap.replaceTop(inserted);
                                maxDistSqr = heap.peek().distSqr;
                            } else {                // heap insert
                                heap.insert(inserted);
                                if(heap.isFull()) {
                                    maxDistSqr = heap.peek().distSqr;
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
                    const Float newOffset = (splitPosition - position[splitAxis]);

                    if(position[splitAxis] < splitPosition) {  //left is closer
                        OpenNode right = actual;
                        actual = OpenNode(actualNode.getLeftChild(), actual.distSqr, actual.distancesSqr);
                        actualNode = nodes[actual.node];
                        right.distSqr -= right.distancesSqr[splitAxis];
                        right.distancesSqr[splitAxis] = newOffset*newOffset;
                        right.distSqr += right.distancesSqr[splitAxis];
                        if(right.distSqr < EPS_FACTOR*maxDistSqr) {
                            right.node = nodes[right.node].getLeftChild() + 1;
                            openNodes.insert(right);
                        }
                    } else {                                 //right is closer
                        OpenNode left = actual;
                        actual = OpenNode(actualNode.getLeftChild()+1, actual.distSqr, actual.distancesSqr);
                        actualNode = nodes[actual.node];
                        left.distSqr -= left.distancesSqr[splitAxis];
                        left.distancesSqr[splitAxis] = newOffset*newOffset;
                        left.distSqr += left.distancesSqr[splitAxis];
                        if(left.distSqr < EPS_FACTOR*maxDistSqr) {
                            left.node = nodes[left.node].getLeftChild();
                            openNodes.insert(left);
                        }
                    }
                }
            }

            return heap.size();
        }
    };
}

#endif