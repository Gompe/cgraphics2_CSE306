#pragma once

#ifndef KD_TREE_H
#define KD_TREE_H

#include <algorithm>
#include <vector>
#include <cassert>

#include "gx_vector.h"
#include "binary_tree.h"

#define KD_TREE_DEBUG 0
#define KD_TREE_CHECK_ALL_LIMIT 10

///////////////////////////////////////////////////////////////////////////////
namespace {
    using namespace binary_tree;
}

template <size_t dim=2>
class KdTree {

    struct SubspaceData {
        size_t first, last;
        size_t axis; 
        double partition_value;

        SubspaceData(const size_t first, const size_t last, const size_t axis, const double partition_value)
            : first(first), last(last), axis(axis), partition_value(partition_value) {}

        SubspaceData(const size_t first, const size_t last)
            : first(first), last(last) {}
    };

    struct KNN {
        // Point P to do the search with
        const Vector P;
        // Number of nearest neighbors to find
        const size_t k;
        // Const Reference to set of points in the KdTree
        const std::vector<Vector> &points;

        // indices of the knn
        std::vector<size_t> knn_indices; 

        KNN(const Vector &P, const size_t k, const std::vector<Vector> &points)
            : P(P), k(k), points(points) {}
        
        double MaxDistance() const {
            if (knn_indices.size() < k) {
                return std::numeric_limits<double>::max();
            }
            else {
                return (points[knn_indices[0]] - P).norm2();
            }
        }

        bool TryInsert(const size_t index) {
            const double distance = (points[index] - P).norm2();
            if (distance >= MaxDistance()) {
                return false;
            }

            auto comp = [this](auto idx1, auto idx2) {
                return (this->P - this->points[idx1]).norm2() < (this->P - this->points[idx2]).norm2();
            };

            // This index should be included in the knn list
            if (knn_indices.size() < k) {
                knn_indices.push_back(index);
            }
            else {
                std::pop_heap(knn_indices.begin(), knn_indices.end(), comp);
                knn_indices[k - 1] = index;
            }

            std::push_heap(knn_indices.begin(), knn_indices.end(), comp);            
            return true;
        } 
    };


    std::vector<Vector> m_points;
    BinaryTree<SubspaceData> m_tree;

public:
    KdTree() = delete;
    KdTree(const KdTree &other) = delete;
    KdTree(KdTree &&other) = delete;
    KdTree& operator=(const KdTree &other) = delete;
    ~KdTree() = default;

    KdTree(const std::vector<Vector> &points)
        : m_points(points) {   

        if (points.size() == 0) {
            m_tree.root = nullptr;
        }
        else {
            m_tree.root = Build(0, points.size(), 0);
        }
    }

    std::vector<Vector> KNearestNeighbors(const Vector &P, const size_t k) const {
        KNN knn_struct(P, k, this->m_points);
        NearestNeighborSearch(m_tree.root, knn_struct);

        std::vector<size_t> indices = knn_struct.knn_indices;
        auto comp = [this, P](auto idx1, auto idx2) {
                return (P - this->m_points[idx1]).norm2() < (P - this->m_points[idx2]).norm2();
            };

        std::sort_heap(indices.begin(), indices.end(), comp);
        std::vector<Vector> out(indices.size());

        std::transform(indices.begin(), indices.end(), out.begin(), [this](auto idx) {
            return this->m_points[idx];
        });

        return out;
    }

    std::vector<Vector> sortedKNearestNeighbors(const Vector &P, int k) const {
        return KNearestNeighbors(P, k);
    }

    void PrintPoints() const {
        std::cout << "--------------------------------------------------------------\n";
        std::cout<< "Printing KdTree - " << m_points.size() << " points\n";
        std::for_each(m_points.begin(), m_points.end(), 
        [](auto p){std::cout << p <<std::endl;});
        std::cout << "--------------------------------------------------------------\n";
    }

private:
    Node<SubspaceData> * Build(const size_t first, const size_t last, const size_t axis) {
        assert(first != last);

        if (last == first + 1) {
            const auto data = SubspaceData(first, last);
            return new Node<SubspaceData>(data);
        }

        // Sort by points by position relative to axis
        auto comp = [axis](auto p1, auto p2) {
            return p1[axis] < p2[axis];
        };

        // Partition the points into two equal halves
        const size_t median = first + (last - first)/2;
        std::nth_element(
            std::next(m_points.begin(), first),
            std::next(m_points.begin(), median),
            std::next(m_points.begin(), last),
            comp
        );
        
        const double partition_value = m_points[median][axis];
        const auto data = SubspaceData(first, last, axis, partition_value);

        const size_t new_axis = (axis + 1) % dim;

        auto node = new Node<SubspaceData>(data);
        node->left = Build(first, median, new_axis);
        node->right = Build(median, last, new_axis);

        return node;
    }

    void NearestNeighborSearch(const Node<SubspaceData> *node, KNN &out) const {
        assert(node != nullptr);
        if (KD_TREE_CHECK_ALL_LIMIT >= node->data.last - node->data.first) {
            // Base case
            for (size_t idx = node->data.first; idx < node->data.last; idx++) {
                out.TryInsert(idx);
            }
        }
        else {
            // Recurse down the tree
            const size_t axis = node->data.axis;
            const double partition_value = node->data.partition_value;

            const Node<SubspaceData> *go_first = node->left, *go_second = node->right;
            if (out.P[axis] > partition_value) {
                std::swap(go_first, go_second);
            }
            
            // Recurse first on this half-space
            NearestNeighborSearch(go_first, out);

            // Recurse on the second branch only if necessary
            const double projlen2 = std::pow(std::abs(out.P[axis] - partition_value), 2.);
            if (projlen2 <= out.MaxDistance()) {
                NearestNeighborSearch(go_second, out);
            }        
        }
    }

};


#endif