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


#ifndef ___POINT_KD_TREE_SHARED__
#define ___POINT_KD_TREE_SHARED__

#include "Config.h"

namespace Importance {

    /// \brief Data structure for storing results of the point kD trees queries. Stores result
    ///        index and its squared distance to the query point
    class KdQueryResult {
    public:
    
        /// \brief index of the result in the array where it is stored
        int index;
    
        union {
        /// \brief squared L2 distance from the query point to the result
        Float distSqr;

        Float weight;
        };        
    
        IMPORTANCE_INLINE KdQueryResult() { }
    
        /// \brief constructs and initialized the query result
        /// \param index index of the result in the array where it is stored
        /// \param distanceSquared squared L2 distance from the query point to the result
        IMPORTANCE_INLINE KdQueryResult(const int index, const Float distanceSquared) : index(index), distSqr(distanceSquared) { }

        IMPORTANCE_INLINE bool operator<(const KdQueryResult& other) const {
            return this->distSqr < other.distSqr;
        }
        IMPORTANCE_INLINE bool operator>(const KdQueryResult& other) const {
            return this->distSqr > other.distSqr;
        }
    };

    class KdQueryResultComparator {
    public:
        IMPORTANCE_INLINE bool operator() (const KdQueryResult& l, const KdQueryResult& r) const {
            return l.distSqr > r.distSqr;   //bude to max-heap
        }
    };

    /// \brief Struct for internal storage of the data. Stores an index of data and its 
    ///        position
    template<class TPoint>
    struct TreeRecordBase {

        /// \brief point position
        TPoint position;

        /// \brief data index
        unsigned int data;

        /// \brief Constructs and initializes an instance
        /// \param data data index
        /// \param position point position
        IMPORTANCE_INLINE TreeRecordBase(const unsigned int data, const TPoint position) : position(position), data(data) { }

        /// \brief Creates default uninitialized instance
        IMPORTANCE_INLINE TreeRecordBase() { }
    };


    /// \brief Opened node during point Kd tree search waiting for processing. Implements the 
    ///        "tracking nodes" (it has saved distances in each axis between its bounding box 
    ///        and the query point) to speed up the search
    template<class TPoint>
    struct OpenKdNodeBase {

        /// \brief Index of the open node
        unsigned int node;

        /// \brief Total squared distance 
        Float distSqr;

        /// \brief Squared distances in each individual axis between the node and query point
        TPoint distancesSqr;

        /// \brief Constructs and initializes the node
        /// \param node Index of the node
        /// \param distSqr Total squared distance 
        /// \param distancesSqr Squared distances in each individual axis between the node 
        ///                     and query point
        IMPORTANCE_INLINE OpenKdNodeBase( const unsigned int node, const Float distSqr, const TPoint distancesSqr ) : 
                node(node), distSqr(distSqr), distancesSqr(distancesSqr) { }

        IMPORTANCE_INLINE OpenKdNodeBase() { }

        IMPORTANCE_INLINE bool operator<(const OpenKdNodeBase& r) const {
            return this->distSqr < r.distSqr;
        }
    };
}

#endif
