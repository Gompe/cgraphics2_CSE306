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

#ifndef VORONOI_NOISE
    #define VORONOI_NOISE 1E-12
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

/// Voronoi Parallel Linear Enumeration
std::vector<Polygon> voronoi(const std::vector<Vector> &voronoiSites){
    double min_x, min_y;
    double max_x, max_y;
    min_x = min_y = std::numeric_limits<double>::max();
    max_x = max_y = std::numeric_limits<double>::lowest();

    for (const auto &P : voronoiSites) {
        min_x = std::min(min_x, P[0]);
        max_x = std::max(max_x, P[0]);

        min_y = std::min(min_y, P[1]);
        max_y = std::max(max_y, P[1]);
    }

    max_x = max_y = std::max(max_x, max_y) + 1;
    min_x = min_y = std::min(min_x, min_y) - 1;

    double numerical_scale = std::abs(max_x) + std::abs(max_y);

    // Defining the enclosing polygon
    Polygon _enclosingPolygon;
    _enclosingPolygon.addVertex(Vector(min_x, min_y, 0.));
    _enclosingPolygon.addVertex(Vector(min_x, max_y, 0.));
    _enclosingPolygon.addVertex(Vector(max_x, max_y, 0.));
    _enclosingPolygon.addVertex(Vector(max_x, min_y, 0.));

    // Make sure that the enclosingPolygon is READ-ONLY
    const Polygon& enclosingPolygon = _enclosingPolygon;

    const int N = voronoiSites.size();
    std::vector<Polygon> voronoiCells(N);

    // #pragma omp parallel for
    for (int siteIdx=0; siteIdx<N; siteIdx++) {
        Vector P = voronoiSites[siteIdx];

        // Add a small random number to P in order to avoid 
        P[0] += VORONOI_NOISE * numerical_scale * random_uniform();
        P[1] += VORONOI_NOISE * numerical_scale * random_uniform();

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

                if ((curVertex - P).norm() < (curVertex - Q).norm()) {
                    if ((prevVertex - P).norm() > (prevVertex - Q).norm()) {
                        outPolygon.addVertex(intersection);
                    }
                    outPolygon.addVertex(curVertex);
                }
                else if ((prevVertex - P).norm() < (prevVertex - Q).norm()) {
                    outPolygon.addVertex(intersection);
                }
            }
            voronoiCell = outPolygon;
        }  

        voronoiCells[siteIdx] = voronoiCell;
    }

    return voronoiCells;
}