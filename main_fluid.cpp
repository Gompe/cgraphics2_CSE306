#include "fluid.h"
#include "gx_polygon.h"
#include "gx_vector.h"
#include "gx_random.h"

int main() {
    std::vector<Vector> points;

    const double box_side = 10;
    Polygon bounding_box;
    bounding_box.AddVertex(Vector(0., 0., 0.));
    bounding_box.AddVertex(Vector(box_side, 0., 0.));
    bounding_box.AddVertex(Vector(box_side, box_side, 0.));
    bounding_box.AddVertex(Vector(0., box_side, 0.));

    Vector box_center = Vector(box_side/2., box_side/2., 0.);

    int N = 5000;
    double x, y;

    for (int i = 0; i < N; i++) {
        boxMuller(1., x, y);
        Vector p = Vector(x, y, 0.) + box_center;

        points.push_back(p);
    }

    Fluid my_fluid = Fluid(points, bounding_box);

    my_fluid.Run(0.05, 100);
}
