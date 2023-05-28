#include "gx_polygon.h"

namespace polygon_detail {
    double MinX(const std::vector<Vector> &points) {
        double out = std::numeric_limits<double>::max();
        for (const auto& point : points) {
            out = std::min(out, point[0]);
        }
        return out;
    }

    double MinY(const std::vector<Vector> &points) {
        double out = std::numeric_limits<double>::max();
        for (const auto& point : points) {
            out = std::min(out, point[1]);
        }
        return out;
    }

    double MaxX(const std::vector<Vector> &points) {
        double out = std::numeric_limits<double>::lowest();
        for (const auto& point : points) {
            out = std::max(out, point[0]);
        }
        return out;
    }

    double MaxY(const std::vector<Vector> &points) {
        double out = std::numeric_limits<double>::lowest();
        for (const auto& point : points) {
            out = std::max(out, point[1]);
        }
        return out;
    }
}

Polygon::Polygon() {}

Polygon::Polygon(const std::vector<Vector> &vertices)
    : vertices(vertices) {}

size_t Polygon::size() const {
    return vertices.size();
}

void Polygon::AddVertex(const Vector &v) {
    vertices.push_back(v);
}

Vector Polygon::operator[](int i) const {
    if (this->size() == 0) {
        fprintf(stderr, "ERROR: Cannot get vertex on empty polygon.\n");
        exit(1);
    }

    i = i % (int) this->size();
    i = (i >= 0) ? i : i + (int) this->size();

    return vertices[i];
} 

Vector& Polygon::operator[](int i) {
    if (this->size() == 0) {
        fprintf(stderr, "ERROR: Cannot get vertex on empty polygon.\n");
        exit(1);
    }

    i = i % (int) this->size();
    i = (i >= 0) ? i : i + (int) this->size();

    return vertices[i];
}

Vector Polygon::GetOutwardsNormal(int i) const {
    const Polygon &self = *this;
    // Normal Vector
    Vector N = Vector(self[i][1] - self[i+1][1], self[i+1][0] - self[i][0], 0.);
    
    return (dot(N, self[i+2] - self[i+1]) < 0) ? N : -N;
}

double Polygon::SignedArea() const {
    double A = 0.;
    const Polygon &self = *this;
    for (int i=0; i < self.size(); i++) {
        A += self[i][0] * self[i+1][1] - self[i+1][0]*self[i][1];
    }
    return 0.5 * A;
}

double Polygon::AbsArea() const {
    return std::abs(SignedArea());
}

Vector Polygon::Centroid() const {
    if (this->size() == 0) {
        fprintf(stderr, "ERROR: Cannot get vertex on empty polygon.\n");
        exit(1);
    }

    double Cx = 0., Cy = 0.;
    const Polygon &self = *this;
    for (int i=0; i < self.size(); i++) {
        double factor = self[i][0] * self[i+1][1] - self[i+1][0]*self[i][1];

        Cx += (self[i][0] + self[i+1][0])*factor;
        Cy += (self[i][1] + self[i+1][1])*factor;
    }

    double factor = 1/(6. * self.SignedArea());
    return Vector(Cx * factor, Cy * factor, 0.);
}

///////////////////////////////////////////////////////////////////////////////
//// Free Functions

double EnergyIntegral(const Polygon &P, const Vector &y) {
    // P has no SignedArea
    if (P.size() < 3) {
        return 0.;
    }

    double integral = 0.;

    std::vector<Vector> c(3);
    c[0] = P[0];

    for (int i = 1; i < P.size() - 1;  i++) {
        c[1] = P[i];
        c[2] = P[i+1];

        // Compute integral inside triangle (c0, c1, c2)
        double sum = 0;
        for (int k = 0; k < 3; k++) {
            for (int l = k; l < 3; l++) {
                sum += dot(c[k] - y, c[l] - y);
            }
        }

        double T = 0.5 * cross(c[1] - c[0], c[2] - c[0]).norm();
        // std::cout << "sum: " << sum << std::endl;

        integral += std::abs(T * sum / 6.);
    }

    if (std::isinf(integral) || std::isnan(integral)) {
        std::cout << "Found weird case:\n";
        std::cout << "Point: " << y << std::endl;
        std::cout << "Polygon: \n";
        for (int i=0; i<P.size(); i++) {
            std::cout << P[i] << "\n";
        }
        std::cout << std::endl;
        exit(1);
    }

    return integral;
}

double OtherEnergyIntegral(const Polygon &P, const Vector &y) {
    // P has no SignedArea
    if (P.size() < 3) {
        return 0.;
    }

    // Pages 102, 103 of lecture notes
    double integral = 0.;
    for (int k=0; k < P.size(); k++) {
        double leftFactor = (P[k-1][0]*P[k][1] - P[k][0]*P[k-1][1]);
        double rightFactor = P[k-1].norm2() + P[k].norm2() + dot(P[k-1], P[k]);

        rightFactor += -4 * dot(y, P[k-1] + P[k]) + 6 * y.norm2();

        integral += leftFactor * rightFactor;
    } 

    integral = std::abs(integral / 12.);

    return integral;
}


Polygon CreateEncompassingSquare(const std::vector<Vector> &points){
    double minCoord = std::min(polygon_detail::MinX(points), polygon_detail::MinY(points));
    double maxCoord = std::max(polygon_detail::MaxX(points), polygon_detail::MaxY(points));

    // Defining the enclosing polygon
    Polygon enclosingSquare;
    enclosingSquare.AddVertex(Vector(minCoord, minCoord, 0.));
    enclosingSquare.AddVertex(Vector(minCoord, maxCoord, 0.));
    enclosingSquare.AddVertex(Vector(maxCoord, maxCoord, 0.));
    enclosingSquare.AddVertex(Vector(maxCoord, minCoord, 0.));

    return enclosingSquare;
}

Polygon CreateDiscretizedDisk(const Vector &x, const double r, const int num_sides) {
    std::vector<Vector> disk_vertices;
    disk_vertices.reserve(num_sides);

    for (int k = 0; k < num_sides; k++) {
        Vector vertex;
        double theta = 2 * M_PI * (double) k / (double) num_sides;

        // cis(theta)
        vertex[0] = cos(theta);
        vertex[1] = sin(theta);

        // scale by radius
        vertex = r * vertex;

        // shift by x
        vertex = x + vertex;

        disk_vertices.push_back(vertex);
    }

    return disk_vertices;
}

std::vector<Polygon> StandardizeCells(
        const std::vector<Polygon> &cells, 
        const Polygon &bounding_box) {
    
    double min_x = polygon_detail::MinX(bounding_box.vertices);
    double min_y = polygon_detail::MinY(bounding_box.vertices);

    double max_x = polygon_detail::MaxX(bounding_box.vertices);
    double max_y = polygon_detail::MaxY(bounding_box.vertices);

    std::vector<Polygon> out_cells;
    out_cells.reserve(cells.size());

    for (const auto &cell : cells) {
        Polygon std_cell;
        for (int i = 0; i < cell.size(); i++) {
            std_cell.AddVertex(
                (cell[i] - Vector(min_x, min_y, 0.)) / Vector((max_x - min_x), (max_y - min_y), 1.)
            );
        }
        out_cells.push_back(std_cell);
    }
    return out_cells;
}

/// Sutherland-Hodgman Polygon Clipping Algorithm
Polygon PolygonClip(Polygon subject, Polygon clip){
    // We pass polygons by value because we need to check that they are in CCW
    // direction and modify them accordingly
    if (subject.SignedArea() < 0) {
        // It is in clockwise direction
        std::reverse(subject.vertices.begin(), subject.vertices.end());
    }
    if (clip.SignedArea() < 0) {
        // It is in clockwise direction
        std::reverse(clip.vertices.begin(), clip.vertices.end());
    }

    // From this point on, we assume that both polygons are in CCW direction
    for (size_t i = 0; i < clip.size(); i++) {
        Polygon result;
        
        // Clip the polygon subject using the edge AB in the clip polygon 
        const Vector &clip_A = clip[i];
        const Vector &clip_B = clip[i+1];

        for (size_t j = 0; j < subject.size(); j++) {
            const Vector &subject_1 = subject[j];
            const Vector &subject_2 = subject[j+1];

            double cross_1 = cross(clip_B - clip_A, subject_1 - clip_A)[2];
            double cross_2 = cross(clip_B - clip_A, subject_2 - clip_A)[2];

            if (cross_1 >= 0 && cross_2 >= 0) {
                // Both points 1 and 2 inside AB edge of clip polygon
                result.AddVertex(subject_2);
            }
            else if (cross_1 >= 0 && cross_2 < 0) {
                // Point 1 inside clip and 2 outside clip
                double t = cross_1 / (cross_1 - cross_2);
                Vector intersection = subject_1 + t * (subject_2 - subject_1);
                result.AddVertex(intersection);
            }
            else if (cross_1 < 0 && cross_2 >= 0) {
                // Point 2 inside clip and 1 outside clip
                double t = cross_1 / (cross_1 - cross_2);
                Vector intersection = subject_1 + t * (subject_2 - subject_1);
                result.AddVertex(intersection);
                result.AddVertex(subject_2);
            }
        }
        
        subject = result;
    }

    return subject;
}

static void DitherVector(Vector &P) {
    P[0] += polygon_detail::dithering_noise * std::abs(P[0]) * (random_uniform() - 0.5);
    P[1] += polygon_detail::dithering_noise * std::abs(P[1]) * (random_uniform() - 0.5);
}

/// Voronoi Parallel Linear Enumeration
std::vector<Polygon> Voronoi(
        const std::vector<Vector> &sites,
        const Polygon &bounding_box) {

    std::vector<double> weights(sites.size());
    std::fill(weights.begin(), weights.end(), 0);
    return FastPowerDiagram(sites, weights, bounding_box);
}

static bool InBoundingBox(const Vector &P, const Polygon &bounding_box) {
    double min_x, min_y, max_x, max_y;
    min_x = polygon_detail::MinX(bounding_box.vertices);
    max_x = polygon_detail::MaxX(bounding_box.vertices);
    min_y = polygon_detail::MinY(bounding_box.vertices);
    max_y = polygon_detail::MaxY(bounding_box.vertices);

    return (P[0] >= min_x - 1E-12) && (P[0] <= max_x + 1E-12) && (P[1] >= min_y - 1E-12) && (P[1] <= max_y + 1E-12);
}

/// Clips the Polygon `subject` using the half-plane containing the points x 
/// such that ||P - x||^2 <= delta_w + ||Q - x||^2
static Polygon ClipByAxis(Polygon subject, const Vector &P, const Vector &Q,
                         double delta_w, const Polygon& bounding_box) {
    Polygon result;
    for (size_t i = 0; i < subject.size(); i++) {
        const Vector &p1 = subject[i];
        const Vector &p2 = subject[i+1];

        bool cond1 = ((p1 - P).norm2() <= delta_w + (p1 - Q).norm2());
        bool cond2 = ((p2 - P).norm2() <= delta_w + (p2 - Q).norm2());

        if (cond1 && cond2) {
            // Both points are inside the half-plane
            result.AddVertex(p2);
        }
        else if (cond1 && !cond2) {
            // Only p1 is inside half-plane
            double t1 = 0.5 * (1 + delta_w/(P-Q).norm2());
            Vector M_prime = P + t1 * (Q - P);

            double t = -dot(p1 - M_prime, P - M_prime)/dot(p2-p1, P - M_prime);
            Vector intersection = p1 + t * (p2 - p1);
            
            if (!InBoundingBox(intersection, bounding_box)) {
                std::cout << "Intersection not in bounding box: " << intersection << std::endl;
                std::cout << "p1 " << p1 << " p2 " << p2 << std::endl;
                std::cout << "P " << P << " Q " << Q << std::endl;
                std::cout << "delta_w " << delta_w << std::endl;
                exit(1);
            }

            result.AddVertex(intersection);
        }
        else if (!cond1 && cond2) {
            // Only p2 is inside half-plane
            double t1 = 0.5 * (1 + delta_w/(P-Q).norm2());
            Vector M_prime = P + t1 * (Q - P);

            double t = -dot(p1 - M_prime, P - M_prime)/dot(p2-p1, P - M_prime);
            Vector intersection = p1 + t * (p2 - p1);
            
            result.AddVertex(intersection);
            result.AddVertex(p2);
        }
    }
    return result;
}

/// power diagram 
std::vector<Polygon> PowerDiagram(
        const std::vector<Vector> &sites,
        const std::vector<double> &weights, 
        const Polygon &bounding_box) {

    const int N = sites.size();
    std::vector<Polygon> cells(N);

    #pragma omp parallel for
    for (int site_idx = 0; site_idx < N; site_idx++) {
        Vector P = sites[site_idx];
        double weight_P = weights[site_idx];

        // Dither P to avoid perfect parallelism and other problematic conditions
        DitherVector(P);

        Polygon cell = bounding_box;

        // Pseudo-Sutherland-Hodgman
        for (int other_idx=0; other_idx < N; other_idx++){
            if (other_idx == site_idx) {
                continue;
            }
            if (cell.size() == 0) {
                break; // Cell won't change anymore
            }
            
            const Vector &Q = sites[other_idx];
            double weight_Q = weights[other_idx];

            double delta_W = weight_P - weight_Q;
            cell = ClipByAxis(cell, P, Q, delta_W, bounding_box);
        }  
        // Store computed result
        cells[site_idx] = cell;
    }
    return cells;
}


std::vector<Polygon> FastPowerDiagram(
        const std::vector<Vector> &sites, 
        const std::vector<double> &weights,
        const Polygon &bounding_box) {

    AcceleratedPowerDiagram accelerationStructure(sites, weights);
    return accelerationStructure.MakeDiagram(bounding_box);
}

///////////////////////////////////////////////////////////////////////////////
//// Accelerated Power Diagram

AcceleratedPowerDiagram::AcceleratedPowerDiagram(
    const std::vector<Vector> &sites,
    const std::vector<double> &weights
    ) : sites(sites), weights(weights), N(sites.size()),
    M(*std::max_element(weights.begin(), weights.end()) + 1.)
{
    std::vector<Vector> liftedSites(N);
    for (int i=0; i<N; i++) {
        liftedSites[i] = sites[i] + Vector(0., 0., sqrt(this->M - weights[i]));
    }
    // Initialize KdTree
    ptrKdTree = std::make_unique< KdTree<3> >(liftedSites);
}

std::vector<Vector> AcceleratedPowerDiagram::LiftedKNN(Vector P, double weight, int k) const {
    P = P + Vector(0., 0., sqrt(M - weight));
    return ptrKdTree->sortedKNearestNeighbors(P, k);
}

std::vector<Polygon> AcceleratedPowerDiagram::MakeDiagram(const Polygon &bounding_box) const {

    std::vector<Polygon> cells(N);

    #pragma omp parallel for
    for (int siteIdx=0; siteIdx<N; siteIdx++) {
        Vector P = sites[siteIdx];
        double weight_P = weights[siteIdx];

        DitherVector(P);
        Vector lifted_P = P + Vector(0., 0., sqrt(M - weight_P));

        Polygon cell = bounding_box;

        int k = polygon_detail::knn_min_k;
        int start_idx = 1;

        std::vector<Vector> nearest_neighbors_3d;

        do {
            // Get k+1 nearest neighbors (nearest neighbor is P itself)
            nearest_neighbors_3d = LiftedKNN(P, weight_P, k+1);

            for (int idx=start_idx; idx < nearest_neighbors_3d.size(); idx++) {
                // Unlift Q
                Vector Q = nearest_neighbors_3d[idx];
                double weight_Q = M - Q[2] * Q[2];
                Q[2] = 0.;

                double delta_W = weight_P - weight_Q;
                cell = ClipByAxis(cell, P, Q, delta_W, bounding_box);

                if (cell.size() == 0) break;
            }

            // Checks to verify if the cell is finished
            if (nearest_neighbors_3d.size() == N) break;
            if (cell.size() == 0) break;

            // Check if it is time to stop because no points can be closer to the cell

            // Distance between lifted P and furthest point in the cell squared
            double largest_cell_dist2 = std::numeric_limits<double>::lowest();
            for (int i=0; i < cell.size(); i++) {
                largest_cell_dist2 = std::min(largest_cell_dist2, (cell[i] - lifted_P).norm2());
            }

            // Distance between P and kth nearest neighbor squared
            double largest_point_dist2 = (lifted_P - nearest_neighbors_3d[nearest_neighbors_3d.size() - 1]).norm2();

            // It is impossible that any point from this point on will change the cell
            if (4 * largest_cell_dist2 <= largest_point_dist2) break;

            // Set new parameters for the search
            start_idx = nearest_neighbors_3d.size();
            k *= 2;
        } while(1);

        // Store computed result
        cells[siteIdx] = cell;
    }
    return cells;
}
