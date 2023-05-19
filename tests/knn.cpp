#include <iostream>
#include "../kd_trees.h"

#include "../gx_random.h"

#include <chrono>

const int N = 1000 * 1000;
const double M = 1000.;

template <typename T>
static void print_vec(const std::vector<T> &vec) {
    for (auto it=vec.begin(); it != vec.end(); it++) {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}

Vector naive_NN(const Vector &P, const std::vector<Vector> &points)
{
    double dist = std::numeric_limits<double>::max();
    Vector out;
    for (const auto &Q : points) {
        if ((P-Q).norm2() < dist) {
            out = Q;
            dist = (P-Q).norm2();
        }
    }
    return out;
}

int main() {
    std::vector<Vector> points;
    for(int i=0; i<N; i++){
        points.push_back(Vector(random_uniform() * M, random_uniform() * M, 0.));
    }

    // 2 dimensional KdTree
    KdTree<2> myTree{points};
    
    const Vector P(random_uniform() * M, random_uniform() * M, 0.);

    std::cout << "Point P: " << P << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    Vector Q = naive_NN(P, points);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Naive - time: " << (end-start).count() << std::endl;
    std::cout << Q << " distance: " << (P-Q).norm() << std::endl;


    start = std::chrono::high_resolution_clock::now();

    Q = myTree.kNearestNeihbors(P, 1)[0];

    end = std::chrono::high_resolution_clock::now();

    std::cout << "Kd Tree - time: " << (end-start).count() << std::endl;

    std::cout << Q << " distance: " << (P-Q).norm() << std::endl;
}