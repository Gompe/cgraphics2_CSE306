#include <iostream>
#include "gx_random.h"
#include "gx_vector.h"
#include "gx_polygon.h"

#include "svg_handler.h"

#include "optimal_transport.h"

#include <random>

int main() {
    std::vector<Vector> points;
    std::vector<double> weights;

    const double eps = 1/16.; 
    for (int i=0; i<20; i++) {
        for(int j=0; j<20; j++){
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

    std::vector<Polygon> cells = powerDiagram(points, weights);

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