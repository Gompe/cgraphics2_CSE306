#include <iostream>
#include "gx_random.h"
#include "gx_vector.h"
#include "gx_polygon.h"

#include "svg_handler.h"

#include <random>

int main() {
    std::vector<Vector> points;


    const double eps = 0; 
    for (int i=0; i<50; i++) {
        for(int j=0; j<50; j++){
            Vector v = Vector(i, j, 0);

            v[0] += random_uniform() * eps;
            v[1] += random_uniform() * eps;

            points.push_back(v);
        }
    }

    // points.push_back(Vector(-1,0,0));
    // points.push_back(Vector(1,0,0));
    // points.push_back(Vector(0, 1, 0));


    std::vector<Polygon> voronoiCells = voronoi(points);

    save_svg(voronoiCells, "out.svg", "blue");
    exit(0);
}