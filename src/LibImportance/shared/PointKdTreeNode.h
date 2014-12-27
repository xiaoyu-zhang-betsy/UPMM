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


#ifndef ___POINT_KD_TREE_NODE__
#define ___POINT_KD_TREE_NODE__

#include "Config.h"

namespace Importance {

    /// \brief Single kD tree node, can be both leaf or inner node. Variant with implicit 
    ///        left child. It  is using packing, so single node has only 8 bytes. 
    ///        
    /// Has to be used in a way that each node's left child is physically stored immediately 
    /// after the parent. When using this type of nodes, the tree can hold the maximum number of 
    /// 2^25 - 1 leaves and cannot have leaves bigger than 127 entries
    class PointKdTreeNode  {
    protected:

        /// \brief First 4 bytes of compact data storage
        ///
        /// Data layout: the int has 32 bits
        /// First 2 bits: type of node/split axis. There are 4 possible values: 
        ///     11 leaf 
        ///     00 X-split inner node
        ///     01 Y-split inner node
        ///     10 Z-split inner node
        /// 
        /// Rest (bits 3 - 32):
        ///     If node is leaf:  first entry index
        ///     If node is inner: index of right child node. Left child is implicitly placed 
        ///                       immediately after this node.
        unsigned int data1;

        /// \brief second 4 bytes of compact data storage.
        /// 
        /// If node is leaf:  bits 1 - 7  size of node
        ///                   bits 8 - 32 index of this node bounding box
        /// If node is inner: split position as Floating point number (needs casting via void* 
        ///                   trick)
        unsigned int data2;

        /// \brief masks anything but highest 2 bits
        static const unsigned int TYPE_MASK = 0xc0000000u;

        /// \brief masks 25 lowest bits
        static const unsigned int LEAF_MASK = 0x01FFFFFFu;

    public:

        static const int MAXIMUM_LEAF_SIZE = 127;

        IMPORTANCE_INLINE PointKdTreeNode () :data1(unsigned(-1)), data2(unsigned(-1)) { }

        /// \brief true if this is leaf node
        /// \return true if this is leaf node
        IMPORTANCE_INLINE bool isLeaf() const {
            return (this->data1 & TYPE_MASK ) == TYPE_MASK;
        }

        /// \brief marks this node as a leaf node
        IMPORTANCE_INLINE void setIsLeaf(){
            this->data1 |= TYPE_MASK;
        }
    
        /// \brief Returns index of first data entry in this leaf. The node must be leaf
        /// \brief return first data entry in this leaf
        IMPORTANCE_INLINE unsigned int getLeafStart() const {
            IMPORTANCE_ASSERT(isLeaf());
            return (this->data1 & ~TYPE_MASK);
        }

        /// \brief sets index of first data entry in this leaf. The node must be leaf
        /// \param index index of first data entry (maximal permitted value is 2^30 - 1)
        IMPORTANCE_INLINE void setLeafStart(const unsigned int index) {
            IMPORTANCE_ASSERT(isLeaf());
            this->data1 &= TYPE_MASK;
            this->data1 |= (index & ~TYPE_MASK);
            IMPORTANCE_ASSERT(isLeaf());
        }

        /// \brief returns number of data entries in the leaf. The node must be leaf
        /// \return number of data entries in the leaf
        IMPORTANCE_INLINE unsigned int getLeafSize() const {
            IMPORTANCE_ASSERT(isLeaf() && ((data2 >> 25) << 25) == (data2&~LEAF_MASK));
            return data2 >> 25;
        }

        /// \brief Sets the number of data entries in this leaf. The node must be leaf
        /// \param size number of data entries in this leaf
        IMPORTANCE_INLINE void setLeafSize(const unsigned int size) {
            IMPORTANCE_ASSERT(isLeaf());
            this->data2 &= LEAF_MASK;
            this->data2 |= (size << 25);
            IMPORTANCE_ASSERT(getLeafSize() == size);
        }

        /// \brief Returns global index of this leaf (used for indexing its bounding box in 
        ///        separate array). The node must be leaf.
        /// \return global index of this leaf
        IMPORTANCE_INLINE unsigned int getLeafNumber() const {
            IMPORTANCE_ASSERT(isLeaf());
            return data2 & LEAF_MASK;
        }

        /// \brief Sets global index of this leaf (used for indexing its bounding box in 
        ///        separate array). The node must be leaf
        /// \param index global index of this leaf (maximum permitted value is 2^25-1)
        IMPORTANCE_INLINE void setLeafNumber(const unsigned int index) {
            IMPORTANCE_ASSERT(isLeaf());
            this->data2 &= (~LEAF_MASK);
            this->data2 |= index;
            IMPORTANCE_ASSERT(getLeafNumber() == index);
        }

        /// \brief Returns axis in which is the node split. Node must be inner
        /// \return axis in which is the node split
        IMPORTANCE_INLINE int getSplitAxis() const {
            IMPORTANCE_ASSERT(!isLeaf() && (this->data1 >> 30) <= 3);
            return this->data1 >> 30;
        }

        /// \brief Returns index of right child node. Node must be inner
        /// \return index of right child node
        IMPORTANCE_INLINE unsigned int getRightChild() const {
            IMPORTANCE_ASSERT(!isLeaf());
            return data1 & ~TYPE_MASK;
        }

        IMPORTANCE_INLINE bool hasRightChild() const {
            return getRightChild() != ~TYPE_MASK;
        }

        /// \brief Sets split axis for this node. Also automatically sets this node as inner node
        /// \brief split axis for this node
        IMPORTANCE_INLINE void setSplitAxis(const unsigned int axis) {
            IMPORTANCE_ASSERT(axis <= 2);
            this->data1 &= ~TYPE_MASK;
            this->data1 |= (axis << 30);
            IMPORTANCE_ASSERT(!isLeaf());
        }

        /// \brief Sets index of right child node. This node must be inner
        /// \param index index of right child (maximal permitted value is 2^30 - 1)
        IMPORTANCE_INLINE void setRightChild(const unsigned int index ) {
            IMPORTANCE_ASSERT(!isLeaf());
            this->data1 &= TYPE_MASK;
            this->data1 |= (index & ~TYPE_MASK);
            IMPORTANCE_ASSERT(!isLeaf());
        }

        /// \brief Returns position of this node's split in world space. Node must be inner
        /// \return position of this node's split
        IMPORTANCE_INLINE float getSplitPosition() const {
            IMPORTANCE_ASSERT(!isLeaf());
            return *((float*)&data2);
        }

        /// \brief Sets position of this node's split in world space. Node must be inner
        /// \param position position of this node's split
        IMPORTANCE_INLINE void setSplitPosition(const float position) {
            IMPORTANCE_ASSERT(!isLeaf());
            this->data2 = *((unsigned int*)&position);
        }
    };


    class Dynamic3dTreeNode {
    protected:

        /// \brief First 4 bytes of compact data storage
        ///
        /// Data layout: the int has 32 bits
        /// First 2 bits: type of node/split axis. There are 4 possible values: 
        ///     11 leaf 
        ///     00 X-split inner node
        ///     01 Y-split inner node
        ///     10 Z-split inner node
        /// 
        /// Rest (bits 3 - 32):
        ///     If node is leaf:  unused
        ///     If node is inner: index of left child node. Right child is implicitly placed 
        ///                       immediately after left child
        unsigned int data1;

        /// \brief second 4 bytes of compact data storage.
        /// 
        /// If node is leaf:  Index of the leaf node array with records
        /// If node is inner: split position as Floating point number (needs casting via void* 
        ///                   trick)
        unsigned int data2;

        /// \brief masks anything but highest 2 bits
        static const unsigned int TYPE_MASK = 0xc0000000u;

    public:

        IMPORTANCE_INLINE Dynamic3dTreeNode () {
            data1 = unsigned int(-1);
            data2 = unsigned int(-1);
        }

        /// \brief true if this is leaf node
        /// \return true if this is leaf node
        IMPORTANCE_INLINE bool isLeaf() const {
            return (this->data1 & TYPE_MASK ) == TYPE_MASK;
        }

        /// \brief marks this node as a leaf node
        IMPORTANCE_INLINE void setIsLeaf(){
            this->data1 |= TYPE_MASK;
        }

        /// \brief Returns index of the array where is the data associated with this leaf node 
        ///        stored. The node must be leaf
        /// \brief return index of the array where is the leaf data stored
        IMPORTANCE_INLINE int getLeafData() const {
            IMPORTANCE_ASSERT(isLeaf());
            return (this->data1 & ~TYPE_MASK);
        }

        /// \brief Sets the index of the array where is the data associated with this leaf node 
        ///        stored. The node must be leaf
        /// \param index index of the array where is the leaf data stored
        IMPORTANCE_INLINE void setLeafData(const int index) {
            IMPORTANCE_ASSERT(isLeaf());
            this->data1 &= TYPE_MASK;
            this->data1 |= (index & ~TYPE_MASK);
            IMPORTANCE_ASSERT(isLeaf());
        }

        /// \brief Returns axis in which is the node split. Node must be inner
        /// \return axis in which is the node split
        IMPORTANCE_INLINE int getSplitAxis() const {
            IMPORTANCE_ASSERT(!isLeaf() && (this->data1 >> 30) <= 3);
            return this->data1 >> 30;
        }

        /// \brief Returns index of left child node. Node must be inner
        /// \return index of left child node
        IMPORTANCE_INLINE int getLeftChild() const {
            IMPORTANCE_ASSERT(!isLeaf());
            return data1 & ~TYPE_MASK;
        }

        /// \brief Returns index of right child node. Node must be inner
        /// \return index of right child node
        IMPORTANCE_INLINE int getRightChild() const {
            return getLeftChild()+1;
        }

        /// \brief Sets index of left child node. This node must be inner
        /// \param index index of left child
        IMPORTANCE_INLINE void setLeftChild(const int index ) {
            IMPORTANCE_ASSERT(!isLeaf());
            this->data1 &= TYPE_MASK;
            this->data1 |= (index & ~TYPE_MASK);
            IMPORTANCE_ASSERT(!isLeaf());
        }

        /// \brief Sets split axis for this node. Also automatically sets this node as inner node
        /// \brief split axis for this node
        IMPORTANCE_INLINE void setSplitAxis(const int axis) {
            IMPORTANCE_ASSERT(axis <= 2);
            this->data1 &= ~TYPE_MASK;
            this->data1 |= (axis << 30);
            IMPORTANCE_ASSERT(!isLeaf());
        }

        /// \brief Returns position of this node's split in world space. Node must be inner
        /// \return position of this node's split
        IMPORTANCE_INLINE float getSplitPosition() const {
            IMPORTANCE_ASSERT(!isLeaf());
            return *((float*)&data2);
        }

        /// \brief Sets position of this node's split in world space. Node must be inner
        /// \param position position of this node's split
        IMPORTANCE_INLINE void setSplitPosition(const float position) {
            IMPORTANCE_ASSERT(!isLeaf());
            this->data2 = *((unsigned int*)&position);
        }
    };
}

#endif