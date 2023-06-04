#include <iostream>
#include <random>

#include <chrono>

#include "gx_polygon.h"
#include "svg_handler.h"

int main() {
    const double box_side = 10;
    Polygon bounding_box;
    bounding_box.AddVertex(Vector(0., 0., 0.));
    bounding_box.AddVertex(Vector(box_side, 0., 0.));
    bounding_box.AddVertex(Vector(box_side, box_side, 0.));
    bounding_box.AddVertex(Vector(0., box_side, 0.));

    Vector box_center = Vector(box_side/2., box_side/2., 0.);

    int N = 10 * 1000;
    double x, y;

    std::vector<Vector> points;
    for (int i = 0; i < N; i++) {
        boxMuller(2., x, y);
        Vector p = Vector(x, y, 0.) + box_center;
        points.push_back(p.clip(0., box_side));
    }

    auto start = std::chrono::high_resolution_clock::now();
    // Actual Code
    auto cells = Voronoi(points, bounding_box);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << "Time elapsed " << diff.count() << "s\n";

    
    for (const auto &pt : points) {
        cells.push_back(CreateDiscretizedDisk(pt, 0.1, 10));
    }

    cells = StandardizeCells(cells, bounding_box);

    // save_svg(cells, "voronoi.svg", "blue");
    exit(0);
}