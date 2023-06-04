#include <iostream>
#include "../kd_trees.h"

#include "../gx_random.h"

#include <chrono>

const int N = 100 * 1000;
const int k = 10;

const double M = 1000.;

template <typename T>
static void print_vec(const std::vector<T> &vec) {
    for (auto it=vec.begin(); it != vec.end(); it++) {
        std::cout << *it << "\n";
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

std::vector<Vector> naive_NN(const Vector &P, const std::vector<Vector> &points, const size_t k) {
    const size_t num_neighbors = std::min(points.size(), k);
    std::vector<std::pair<double, Vector>> out(num_neighbors);

    std::fill(out.begin(), out.end(), std::make_pair(std::numeric_limits<double>::max(), Vector()));

    for (const auto &Q : points) {
        const double dist = out[0].first;

        if ((P-Q).norm2() < dist) {
            std::pop_heap(out.begin(), out.end(), [](auto a, auto b){
                return a.first < b.first;
            });

            out[num_neighbors - 1] = std::make_pair((P-Q).norm2(), Q);
            std::push_heap(out.begin(), out.end(), [](auto a, auto b){
                return a.first < b.first;
            });
        }
    }

    std::vector<Vector> true_out;
    for (auto p : out) {
        true_out.push_back(p.second);
    }

    std::sort_heap(true_out.begin(), true_out.end(), [&P](auto Q, auto R){
        return (P - Q).norm2() < (P - R).norm2();
    });

    return true_out;
}

int main() {
    std::vector<Vector> points;
    for(int i=0; i<N; i++){
        points.push_back(Vector(random_uniform() * M, random_uniform() * M, 0.));
    }

    // 2 dimensional KdTree
    KdTree<3> myTree{points};
    
    const Vector P(random_uniform() * M, random_uniform() * M, 0.);

    std::cout << "Point P: " << P << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto list_naive = naive_NN(P, points, k);

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Naive - time: " << (end-start).count() << std::endl;

    start = std::chrono::high_resolution_clock::now();

    auto list_kd = myTree.sortedKNearestNeighbors(P, k);

    end = std::chrono::high_resolution_clock::now();

    std::cout << "Kd Tree - time: " << (end-start).count() << std::endl;

    auto printer = [&P](auto l) {
        for (auto Q : l) {
            std::cout << Q << " " << (P - Q).norm2() << std::endl;
        }
        std::cout << std::endl;
    };

    printer(list_kd);
    std::cout << "------------------------------------------------------------\n";

    printer(list_naive);

}