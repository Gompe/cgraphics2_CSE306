#include "laguerre.h"
#include "svg_handler.h"

#include <iostream>
#include <random>

int main() {
    const double box_side = 10;
    Polygon bounding_box;
    bounding_box.AddVertex(Vector(0., 0., 0.));
    bounding_box.AddVertex(Vector(box_side, 0., 0.));
    bounding_box.AddVertex(Vector(box_side, box_side, 0.));
    bounding_box.AddVertex(Vector(0., box_side, 0.));

    Vector box_center = Vector(box_side/2., box_side/2., 0.);

    int N = 5;
    double x, y;

    std::vector<Vector> points;
    for (int i = 0; i < N; i++) {
        boxMuller(1., x, y);
        Vector p = Vector(x, y, 0.) + box_center;
        points.push_back(p.clip(0., box_side));
    }


    std::vector<double> weights(N);

    auto cells = FastPowerDiagram(points, weights, bounding_box);

    std::vector<double> lambdas(cells.size());
    std::fill(lambdas.begin(), lambdas.end(), bounding_box.AbsArea()/(double) cells.size());

    std::vector<Polygon> out;
    out.insert(out.begin(), cells.begin(), cells.end());
    for (auto point : points) {
        out.push_back(CreateDiscretizedDisk(point, 0.05, 50));
    }
    out = StandardizeCells(out, bounding_box);
    save_svg(out, "laguerre_before.svg", "blue");

    // Create Laguerre Structure
    LaguerreDiagram my_diagram(points, cells, lambdas, weights, bounding_box);
    my_diagram.Optimize();

    std::cout << "Weights\n";
    for (int i = 0; i < my_diagram.cells.size(); i++) {
        std::cout << my_diagram.weights[i] << " ";
    }
    std::cout << "\n";

    out.clear();
    out.insert(out.begin(), my_diagram.cells.begin(), my_diagram.cells.end());
    for (auto point : points) {
        out.push_back(CreateDiscretizedDisk(point, 0.05, 50));
    }
    out = StandardizeCells(out, bounding_box);

    save_svg(out, "laguerre_after.svg", "blue");

    exit(0);
}