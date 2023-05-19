#include <iostream>
#include <chrono>

#include "gx_random.h"
#include "gx_vector.h"
#include "gx_polygon.h"

#include "svg_handler.h"

#include "optimal_transport.h"

#include "kd_trees.h"

#include <random>

int main() {
    std::vector<Vector> points;
    std::vector<double> weights;

    int N = 100;

    const double eps = 1/16.; 
    for (int i=0; i<N; i++) {
        for(int j=0; j<N; j++){
            Vector v = Vector(i * random_uniform(), j * random_uniform(), 0);
            weights.push_back(0);
            v[0] += random_uniform() * eps;
            v[1] += random_uniform() * eps;

            points.push_back(v);
        }
    }

    // points.push_back(Vector(-1,0,0));
    // points.push_back(Vector(1,0,0));
    // points.push_back(Vector(0, 1, 0));  

    std::cout << "naive power diagram" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Polygon> cells = powerDiagram(points, weights);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "time: " << (end - start).count()/(double) 1E6 << std::endl;
    std::cout << std::endl;

    std::cout << "fast power diagram" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    cells = fastPowerDiagram(points, weights);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "time: " << (end - start).count()/(double) 1E6 << std::endl;

    std::vector<double> lambdas(cells.size());
    std::fill(lambdas.begin(), lambdas.end(), 1/(double) cells.size());

    save_svg(cells, "before.svg", "blue");

    optimalTransport(
        points,
        cells,
        lambdas,
        weights
    );

    save_svg(cells, "after.svg", "blue");

    exit(0);
}