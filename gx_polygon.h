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

/// Voronoi Parallel Linear Enumeration
std::vector<Polygon> voronoi(const std::vector<Vector> &voronoiSites){

    const Polygon enclosingPolygon = createEncompassingSquare(voronoiSites);
    double dithering_scale = std::abs(enclosingPolygon[0][0]) + std::abs(enclosingPolygon[1][0]);

    const int N = voronoiSites.size();
    std::vector<Polygon> voronoiCells(N);

    #pragma omp parallel for
    for (int siteIdx=0; siteIdx<N; siteIdx++) {
        Vector P = voronoiSites[siteIdx];

        // Add a small random number to P in order to avoid 
        P[0] += DITHERING_NOISE * dithering_scale * random_uniform();
        P[1] += DITHERING_NOISE * dithering_scale * random_uniform();

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
        double weightP = weights[siteIdx];

        // Add a small random number to P in order to avoid 
        P[0] += DITHERING_NOISE * dithering_scale * random_uniform();
        P[1] += DITHERING_NOISE * dithering_scale * random_uniform();

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