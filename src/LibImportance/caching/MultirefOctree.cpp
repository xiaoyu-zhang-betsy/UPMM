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


#include "MultirefOctree.h"
#include "../shared/StaticArray.h"
#include "MultirefOctreeAux.h"
using namespace Importance;

void MultirefOctree::iterate( const int nodeIndex, const BoundingBox3& bBox, const int depth ) {
    const Vector3f center = bBox.getCenter();
    nodes[nodeIndex].center = center;
    
    if(nodes[nodeIndex].records.size() <= this->leafSize || depth >= this->maxDepth) {
        nodes[nodeIndex].setLeaf(true);
        return;
    }

    int childOffset = int(this->nodes.size());
    nodes[nodeIndex].setLeaf(false);
    nodes[nodeIndex].setChildPtr(childOffset);
    for(int i = 0; i < 8; ++i) {
        this->nodes.push(Node());
    }

    // vector might get reallocated
    Node& actual2 = this->nodes[nodeIndex];
    int ownPtr = 0;
    ImStaticArray<BoundingBox3, 8> childBboxes;
    this->createChildBboxes(bBox, childBboxes);
    
    for(int i = 0; i < actual2.records.size(); ++i) {
        const float validity = this->records[actual2.records[i]].radius();
        const Vector3f recordPos = this->records[actual2.records[i]].position();
        BoundingSphere sphere(recordPos, validity);
        if(nonSplittable(sphere,bBox)) {  // leave it in this node
            actual2.records[ownPtr++] = actual2.records[i];
        } else {
            for(int j = 0; j < 8; ++j) {
                if(sphere.intersect(childBboxes[j])) {
                        this->nodes[childOffset+j].records.push(actual2.records[i]);
                }
            }
        }
    }

    actual2.records.resize(ownPtr);

    for(int i = 0; i < 8; ++i) {
        iterate(childOffset+i, childBboxes[i], depth+1);
    }
}

void MultirefOctree::addRecord(const Vector3f position, const float radius, const float validityRadiusMult) {
    this->records.push(Record(position, /*normal,*/ radius*validityRadiusMult));
    this->updateMaxDepth();
    addRecordRecursively(0, int(this->records.size()-1), entireBbox, 0);
}

void MultirefOctree::addRecordRecursively( const int nodeIndex, const int recordId, 
                                             const BoundingBox3& bBox, const int depth ) {
    const Record& record = this->records[recordId];
    Node& node = this->nodes[nodeIndex];

    if(node.isLeaf()) {
        node.records.push(recordId);
        if(node.records.size() > this->leafSize && depth < this->maxDepth) {
            //split
            Vector3f center = bBox.getCenter();
            
            int childOffset = int(this->nodes.size());
            node.setLeaf(false);
            node.setChildPtr(childOffset);
            node.center = center;
            ImStaticArray<BoundingBox3, 8> childBboxes;
            this->createChildBboxes(bBox, childBboxes);

            for(int i = 0; i < 8; ++i) {
                this->nodes.push(Node());
                this->nodes.back().setLeaf(true);
                this->nodes.back().records.reserve(this->leafSize);
            }

            // vector might get reallocated
            int ownPtr = 0;
            Node& node2 = this->nodes[nodeIndex];

            for(int i = 0; i < node2.records.size(); ++i) {
                const float validity = this->records[node2.records[i]].radius();
                const Vector3f recordPos = this->records[node2.records[i]].position();
                BoundingSphere sphere(recordPos, validity);
                if(nonSplittable(sphere, bBox)) {   // leave it in this node
                    node2.records[ownPtr++] = node2.records[i];
                } else {
                    for(int j = 0; j < 8; ++j) {
                        if(sphere.intersect(childBboxes[j])) {
                            this->nodes[childOffset+j].records.push(node2.records[i]);
                        }
                    }
                }
            }
            node2.records.resize(ownPtr);
        }
    } else {
        const int children = node.getChildPtr();
        const BoundingSphere sphere(record.position(), record.radius());
        if(nonSplittable(sphere, bBox)) {
            node.records.push(recordId);
        } else {

            ImStaticArray<BoundingBox3, 8> childBboxes;
            this->createChildBboxes(bBox, childBboxes);

            for(int i = 0; i < 8; ++i) {
                if(sphere.intersect(childBboxes[i])) {
                    addRecordRecursively(children+i, recordId, childBboxes[i], depth+1);
                }
            }
        }
    }
}

//void MultirefOctree::makeStats() const {
//    int totalReferences = 0;
//    //int nonLeafReferences = 0;
//
//    for(int i = 0; i < this->nodes.size(); ++i) {
//        totalReferences += this->nodes[i].records.size();
//
//        //if(!this->nodes[i].isLeaf()) {
//        //    nonLeafReferences += this->nodes[i].records.size();
//        //}
//    }
//
//    std::ofstream os("C:/multirefOctree.txt");
//    std::map<int, int> nodes;
//    std::map<int, int> records;
//    statsIterate(nodes, records, 0, 0);
//    os << "Nodes in each level:" << std::endl;
//    for(auto it = nodes.begin(); it != nodes.end(); ++it) {
//        os << it->first << ": " << it->second << std::endl;
//    }
//    os << std::endl << std::endl;
//
//    os << "Records in each level:" << std::endl;
//    for(auto it = records.begin(); it != records.end(); ++it) {
//        os << it->first << ": " << it->second << std::endl;
//    }
//    os << std::endl << std::endl;
//    os << "Records total: " << this->records.size();
//}