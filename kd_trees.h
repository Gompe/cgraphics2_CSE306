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

using namespace binary_tree;
using RandomIt = std::vector<Vector>::iterator;

/// KdTree for the given number of dimensions
template <size_t dim=2>
class KdTree
{
    // Point Set
    struct PS {
        // Value at which the points are partitioned
        double medianValue;
        // indexes of first and last points in the set
        RandomIt first, last;
        // axis on which the set is partitioned
        int axis;
    };

    struct NearestNeighborStruct {
        // K Nearest Neighbors of P
        const Vector P;
        const int k;
        std::vector<std::pair<double, RandomIt>> knn;

        NearestNeighborStruct(const Vector &P, int k) : P(P), k(k){
            knn.reserve(k);
        }

        double maxDistance() const {
            if (knn.size() < k)
                return std::numeric_limits<double>::max();
            return knn[0].first;
        }

        bool tryInsert(RandomIt it) {
            double distance = (P - *it).norm2();
            if (distance >= maxDistance()) {
                return false;
            }
            
            if (knn.size() < k) {
                knn.push_back(std::make_pair(distance, it));
            }
            else {
                std::pop_heap(knn.begin(), knn.end());
                knn[k-1] = std::make_pair(distance, it);
            }

            std::push_heap(knn.begin(), knn.end());
            return true;
        }
    };

    std::vector<Vector> points;
    BinaryTree<PS> tree;

public:
    KdTree() = delete;
    KdTree(const KdTree &other) = delete;
    KdTree& operator=(const KdTree &other) = delete;

    KdTree(const std::vector<Vector> &points)
    : points(points)
    {   
        tree.root = new Node<PS>();
        build(this->points.begin(), this->points.end(), 0, tree.root);
    }

    /// Returns the k nearest neighbors to point P in the set points
    /// in an ARBITRARY ORDER
    std::vector<Vector> kNearestNeihbors(const Vector &P, int k) const
    {   
        // Initialize structure
        NearestNeighborStruct nns(P, k);

        // Search for the nearest neighbors
        nearestNeighborSearch(tree.root, nns);

        std::vector<Vector> output;
        std::for_each(nns.knn.begin(), nns.knn.end(), 
            [&output](const auto &p){output.push_back(*p.second);});

        return output;
    }

    /// Returns the k nearest neighbors to point P in the set points
    /// sorted by distance.
    /// Ties are decided aribitrarily.
    std::vector<Vector> sortedKNearestNeighbors(const Vector &P, int k) const
    {
        auto knn = kNearestNeihbors(P, k);
        std::stable_sort(knn.begin(), knn.end(), 
            [&P](const Vector &Q, const Vector &R){ return (P-Q).norm2() < (P-R).norm2(); }
            );
        
        return knn;
    }

private:
    void build(RandomIt first, RandomIt last, int axis, Node<PS> *node)
    {

        #if KD_TREE_DEBUG
            assert(first != last);
        #endif

        #if KD_TREE_DEBUG
            printf("Build Trace: first index=%d, last index=%d\n",
                (int) std::distance<std::vector<Vector>::const_iterator>(points.begin(), first),
                (int) std::distance<std::vector<Vector>::const_iterator>(points.begin(), last)
                );
        #endif
        
        node->data.first = first;
        node->data.last = last;
        node->data.axis = axis;

        int n = std::distance(first, last);
        if (n == 1) {
            return;
        }

        RandomIt median = std::next(first, n/2);

        // Compares a < b along axis
        auto comparer = [axis](const Vector &a, const Vector &b){
            return (a[axis] < b[axis]);
        };

        // Partial Sort where median element goes to the right place
        std::nth_element(first, median, last, comparer);
        node->data.medianValue = (*median)[axis];

        int newAxis = (axis + 1) % dim;

        node->left = new Node<PS>();
        build(first, median, newAxis, node->left);

        node->right = new Node<PS>();
        build(median, last, newAxis, node->right);
    }

    void nearestNeighborSearch(const Node<PS> *node, NearestNeighborStruct &nns) const
    {   
        // Base Cases -- there are few nearest neighbors to check
        if (node == nullptr) {
            return; 
        }

        #if KD_TREE_DEBUG
            printf("Search Trace: first index=%d, last index=%d\n",
                (int) std::distance<std::vector<Vector>::const_iterator>(points.begin(), node->data.first),
                (int) std::distance<std::vector<Vector>::const_iterator>(points.begin(), node->data.last)
                );
        #endif

        if (KD_TREE_CHECK_ALL_LIMIT >= std::distance(node->data.first, node->data.last)) {
            for (RandomIt it=node->data.first; it != node->data.last; it++) {
                nns.tryInsert(it);
            }
            return;
        }


        int axis = node->data.axis;
        const Node<PS> *goFirst = node->left, *goSecond = node->right;
        if (nns.P[axis] > node->data.medianValue) 
            std::swap(goFirst, goSecond);

        nearestNeighborSearch(goFirst, nns);
        // Check if we have to check the other branch
        if (std::abs(node->data.medianValue - nns.P[axis]) < nns.maxDistance())
            nearestNeighborSearch(goSecond, nns);
    }
};

#endif