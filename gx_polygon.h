/**
 * gx_polygon.h
 * 
 *  Implements the Polygon logic. Polygons are implemented in 2 dimensions.
 * Even though we use Vector, we limit ourselfs to vectors with z-component
 * equal to 0.
*/

#pragma once

#include "gx_vector.h"
#include "gx_random.h"

// Nearest Neighbors acceleration
#include "kd_trees.h"

#include <memory>
#include <cassert>

#ifndef EPSILON_ZERO
    #define EPSILON_ZERO 1E-9
#endif

#ifndef EPSILON
    #define EPSILON 1E-4
#endif

#ifndef DITHERING_NOISE
    #define DITHERING_NOISE 1E-12
#endif

// Number of nearest neighbors to start with
#ifndef KNN_MIN_K
    #define KNN_MIN_K 20
#endif

class Polygon {  
public:
    Polygon() {}

    Polygon(const std::vector<Vector> &vertices)
    : vertices(vertices) {}

    std::vector<Vector> vertices;

    size_t size() const {
        return vertices.size();
    }

    void addVertex(const Vector &v) {
        vertices.push_back(v);
    }

    /// polygon[i] returns the vertex given by i mod number of vertices
    Vector operator[](int i) const {
        if (this->size() == 0) 
            throw std::logic_error("Cannot get vertex on empty polygon.\n");

        i = i % (int) this->size();
        i = (i >= 0) ? i : i + (int) this->size();

        return vertices[i];
    } 

    /// Returns an outwards normal of the edge connecting vertex[i] to vertex[i+1]
    Vector getOutwardsNormal(int i) const {
        const Polygon &self = *this;
        // Normal Vector
        Vector N = Vector(self[i][1] - self[i+1][1], self[i+1][0] - self[i][0], 0.);
        
        return (dot(N, self[i+2] - self[i+1]) < 0) ? N : -N;
    }

    /// Returns the directed area of the polygon
    double area() const {
        double A = 0.;
        const Polygon &self = *this;
        for (int i=0; i < self.size(); i++) {
            A += self[i][0] * self[i+1][1] - self[i+1][0]*self[i][1];
        }
        return 0.5 * A;
    }

    /// Returns (absolute) area of the polygon
    double absArea() const {
        return std::abs(area());
    }

    /// Returns the centroid of the polygon (Assumes it is not self-intersecting)
    Vector centroid() const {

        if (this->size() == 0)
            throw std::logic_error("ERROR: Empty polygon does not have centroid.\n");

        double Cx = 0., Cy = 0.;
        const Polygon &self = *this;
        for (int i=0; i < self.size(); i++) {
            double factor = self[i][0] * self[i+1][1] - self[i+1][0]*self[i][1];

            Cx += (self[i][0] + self[i+1][0])*factor;
            Cy += (self[i][1] + self[i+1][1])*factor;
        }

        double factor = 1/(6 * self.area());
        return Vector(Cx * factor, Cy * factor, 0.);
    }
};  

/// Returns the real number t for which A + t*(B-A) is inside the line segment CD
static double edgeIntersection(const Vector &A, const Vector &B, const Vector &C, const Vector &D){
    // Outwards normal to CD
    Vector N = Vector(C[1] - D[1], D[0] - C[0], 0.);
    
    double denominator = dot(B-A, N);
    if (std::abs(denominator) < EPSILON_ZERO)
        return std::numeric_limits<double>::max();
    
    double numerator = dot(C - A, N);
    return numerator/denominator;
}

/// Verifies if a point P is inside the halfplane defined by line passing 
/// through A with outward normal N
static bool testInside(const Vector &P, const Vector &A, const Vector &N){
    return (dot(P-A, N) <= EPSILON);
}


/// Sutherland-Hodgman Polygon Clipping Algorithm
Polygon polygonClip(Polygon subjectPolygon, const Polygon &clipPolygon){
    for(int edgeIdx=0; edgeIdx < clipPolygon.size(); edgeIdx++){
        Polygon outPolygon = Polygon();

        Vector prevClipVertex = clipPolygon[edgeIdx];
        Vector currClipVertex = clipPolygon[edgeIdx + 1];
        Vector outwardNormal = clipPolygon.getOutwardsNormal(edgeIdx);

        for(int i=0; i<subjectPolygon.size(); i++){
            Vector curVertex = subjectPolygon[i];
            Vector prevVertex = subjectPolygon[i-1];

            double t = edgeIntersection(prevVertex, curVertex, prevClipVertex, currClipVertex);
            Vector intersection = prevVertex + t * (curVertex - prevVertex);

            bool isCurVertexInside = testInside(curVertex, prevClipVertex, outwardNormal);
            bool isPrevVertexInside = testInside(prevVertex, prevClipVertex, outwardNormal);

            if (isCurVertexInside) {
                if (!isPrevVertexInside) {
                    outPolygon.addVertex(intersection);
                }
                outPolygon.addVertex(curVertex);
            } 
            else if (isPrevVertexInside) {
                outPolygon.addVertex(intersection);
            }
        }
        subjectPolygon = outPolygon;
    }

    return subjectPolygon;
}

/// Creates a square encompassing all points in the container of Vectors
Polygon createEncompassingSquare(const std::vector<Vector> &points){
    double minCoord = std::accumulate(points.begin(), points.end(), 
        std::numeric_limits<double>::max(), [](double acc, const Vector P) {
            return std::min(acc, std::min(P[0], P[1]));
        }) - 1.;
    
    double maxCoord = std::accumulate(points.begin(), points.end(), 
        std::numeric_limits<double>::min(), [](double acc, const Vector P) {
            return std::max(acc, std::max(P[0], P[1]));
        }) + 1.;

    // Defining the enclosing polygon
    Polygon enclosingSquare;
    enclosingSquare.addVertex(Vector(minCoord, minCoord, 0.));
    enclosingSquare.addVertex(Vector(minCoord, maxCoord, 0.));
    enclosingSquare.addVertex(Vector(maxCoord, maxCoord, 0.));
    enclosingSquare.addVertex(Vector(maxCoord, minCoord, 0.));

    return enclosingSquare;
}

static void ditherVector(Vector &P) 
{
    P[0] += DITHERING_NOISE * std::abs(P[0]) * (random_uniform() - 0.5);
    P[1] += DITHERING_NOISE * std::abs(P[1]) * (random_uniform() - 0.5);
}

/// Voronoi Parallel Linear Enumeration
std::vector<Polygon> voronoi(const std::vector<Vector> &voronoiSites){

    const Polygon enclosingPolygon = createEncompassingSquare(voronoiSites);

    const int N = voronoiSites.size();
    std::vector<Polygon> voronoiCells(N);

    #pragma omp parallel for
    for (int siteIdx=0; siteIdx<N; siteIdx++) {
        Vector P = voronoiSites[siteIdx];
        // Add small noise to avoid points in general position
        ditherVector(P);

        Polygon voronoiCell = enclosingPolygon;

        // Pseudo-Sutherland-Hodgman
        for (int otherSiteIdx=0; otherSiteIdx < N; otherSiteIdx++){
            if (otherSiteIdx == siteIdx) 
                continue;
            
            const Vector &Q = voronoiSites[otherSiteIdx];
            // Use the bissector of PQ to cut the voronoi cell

            Polygon outPolygon = Polygon();

            Vector M = 0.5*(P+Q);

            for(int i=0; i<voronoiCell.size(); i++){
                const Vector& curVertex = voronoiCell[i];
                const Vector& prevVertex = voronoiCell[i-1];

                double t = dot(M - prevVertex, Q - P)/dot(curVertex - prevVertex, Q - P);
                Vector intersection = prevVertex + t * (curVertex - prevVertex);

                bool condCurVertex = ((curVertex - P).norm2() < (curVertex - Q).norm2());
                bool condPrevVertex = ((prevVertex - P).norm2() < (prevVertex - Q).norm2());

                if (condCurVertex) {
                    if (!condPrevVertex) {
                        outPolygon.addVertex(intersection);
                    }
                    outPolygon.addVertex(curVertex);
                }
                else if (condPrevVertex) {
                    outPolygon.addVertex(intersection);
                }
            }
            voronoiCell = outPolygon;
        }  

        // Store the computed result
        voronoiCells[siteIdx] = voronoiCell;
    }
    return voronoiCells;
}

/// power diagram 
std::vector<Polygon> powerDiagram(const std::vector<Vector> &sites, const std::vector<double> &weights){

    const Polygon enclosingPolygon = createEncompassingSquare(sites);
    double dithering_scale = std::abs(enclosingPolygon[0][0]) + std::abs(enclosingPolygon[1][0]);

    const int N = sites.size();
    std::vector<Polygon> powerDiagramCells(N);

    #pragma omp parallel for
    for (int siteIdx=0; siteIdx<N; siteIdx++) {
        Vector P = sites[siteIdx];
        ditherVector(P);
        double weightP = weights[siteIdx];

        Polygon cell = enclosingPolygon;

        // Pseudo-Sutherland-Hodgman
        for (int otherSiteIdx=0; otherSiteIdx < N; otherSiteIdx++){
            if (otherSiteIdx == siteIdx) 
                continue;
            
            if (cell.size() == 0) {
                // Cell cannot get smaller
                break;
            }

            const Vector &Q = sites[otherSiteIdx];
            double weightQ = weights[otherSiteIdx];

            // Use the bissector of PQ to cut the voronoi cell

            Polygon outPolygon = Polygon();

            Vector M = 0.5*(P+Q);
            Vector Mprime = M + (Q - P) * (weightP - weightQ)/(2*(Q - P).norm2());

            for(int i=0; i<cell.size(); i++){
                const Vector& curVertex = cell[i];
                const Vector& prevVertex = cell[i-1];

                double t = dot(Mprime - prevVertex, Q - P)/dot(curVertex - prevVertex, Q - P);
                Vector intersection = prevVertex + t * (curVertex - prevVertex);

                bool condCurVertex = ((curVertex - P).norm2() - weightP < (curVertex - Q).norm2() - weightQ);
                bool condPrevVertex = ((prevVertex - P).norm2() - weightP < (prevVertex - Q).norm2() - weightQ);

                if (condCurVertex) {
                    if (!condPrevVertex) {
                        outPolygon.addVertex(intersection);
                    }
                    outPolygon.addVertex(curVertex);
                }
                else if (condPrevVertex) {
                    outPolygon.addVertex(intersection);
                }
            }
            cell = outPolygon;
        }  
        // Store computed result
        powerDiagramCells[siteIdx] = cell;
    }
    return powerDiagramCells;
}


class AcceleratedPowerDiagram 
{
    const std::vector<Vector> sites;
    const std::vector<double> weights;
    const int N;
    const double M;

    std::unique_ptr< KdTree<3> > ptrKdTree;

public:
    AcceleratedPowerDiagram& operator=(const AcceleratedPowerDiagram &other) = delete;
    AcceleratedPowerDiagram(const AcceleratedPowerDiagram &other) = delete;

    AcceleratedPowerDiagram(
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

    std::vector<Vector> liftedKNN(Vector P, double weight, int k) const
    {
        P = P + Vector(0., 0., sqrt(M - weight));
        return ptrKdTree->sortedKNearestNeighbors(P, k);
    }

    std::vector<Polygon> makeDiagram() const
    {
        const Polygon enclosingPolygon = createEncompassingSquare(sites);
        double dithering_scale = std::abs(enclosingPolygon[0][0]) + std::abs(enclosingPolygon[1][0]);

        std::vector<Polygon> cells(N);

        // #pragma omp parallel for
        for (int siteIdx=0; siteIdx<N; siteIdx++) {
            Vector P = sites[siteIdx];
            ditherVector(P);

            double weightP = weights[siteIdx];
            Vector liftedP = P + Vector(0., 0., sqrt(M - weightP));

            Polygon cell = enclosingPolygon;

            int k = KNN_MIN_K, searchStart = 1;
            std::vector<Vector> liftedNearestNeighbors;

            do {
                // Get k+1 nearest neighbors (nearest neighbor is P itself)
                liftedNearestNeighbors = liftedKNN(P, weightP, k+1);

                for (int idx=searchStart; idx < liftedNearestNeighbors.size(); idx++) {
                    // Unlift Q
                    Vector Q = liftedNearestNeighbors[idx];
                    double weightQ = M - Q[2] * Q[2];
                    Q[2] = 0.;

                    Vector M = 0.5*(P+Q);
                    Vector Mprime = M + (Q - P) * (weightP - weightQ)/(2*(Q - P).norm2());

                    Polygon outPolygon = Polygon();
                    for(int i=0; i<cell.size(); i++){
                        const Vector& curVertex = cell[i];
                        const Vector& prevVertex = cell[i-1];

                        double t = dot(Mprime - prevVertex, Q - P)/dot(curVertex - prevVertex, Q - P);
                        Vector intersection = prevVertex + t * (curVertex - prevVertex);

                        bool condCurVertex = ((curVertex - P).norm2() - weightP < (curVertex - Q).norm2() - weightQ);
                        bool condPrevVertex = ((prevVertex - P).norm2() - weightP < (prevVertex - Q).norm2() - weightQ);

                        if (condCurVertex) {
                            if (!condPrevVertex) {
                                outPolygon.addVertex(intersection);
                            }
                            outPolygon.addVertex(curVertex);
                        }
                        else if (condPrevVertex) {
                            outPolygon.addVertex(intersection);
                        }
                    }
                    cell = outPolygon;
                    if (cell.size() == 0) break;
                }

                // Checks to verify if the cell is finished
                if (liftedNearestNeighbors.size() == N) break;
                if (cell.size() == 0) break;

                // Check if it is time to stop because no points can be closer to the cell

                // Distance between lifted P and furthest point in the cell squared
                double distanceFurthestPoint2 = std::numeric_limits<double>::lowest();
                for (int i=0; i < cell.size(); i++) {
                    distanceFurthestPoint2 = std::min(distanceFurthestPoint2, (cell[i] - liftedP).norm2());
                }

                // Distance between P and kth nearest neighbor squared
                double distanceFurthestSite2 = (liftedP - liftedNearestNeighbors[liftedNearestNeighbors.size() - 1]).norm2();

                // It is impossible that any point from this point on will change the cell
                if (4 * distanceFurthestPoint2 <= distanceFurthestSite2) break;

                // Set new parameters for the search
                searchStart = liftedNearestNeighbors.size();
                k *= 2;
            } while(1);

            // Store computed result
            cells[siteIdx] = cell;
        }
        return cells;
    }
};

std::vector<Polygon> fastPowerDiagram(const std::vector<Vector> &sites, const std::vector<double> &weights)
{
    AcceleratedPowerDiagram accelerationStructure(sites, weights);
    return accelerationStructure.makeDiagram();
}