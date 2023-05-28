#pragma once

#ifndef POLYGON_H
#define POLYGON_H

#include "gx_vector.h"
#include "gx_random.h"

// Nearest Neighbors acceleration
#include "kd_trees.h"

#include <memory>
#include <cassert>

namespace polygon_detail {
    const double numeric_epsilon = 1E-12;
    const double dithering_noise = 1E-12;

    const int knn_min_k = 20;
    const int num_sides_circle = 50;

    double MinX(const std::vector<Vector> &points);
    double MinY(const std::vector<Vector> &points);
    double MaxX(const std::vector<Vector> &points);
    double MaxY(const std::vector<Vector> &points);
    
}; // namespace polygon_detail

class Polygon {  
public:
    std::vector<Vector> vertices;

    Polygon();
    Polygon(const std::vector<Vector> &vertices);

    size_t size() const;

    void AddVertex(const Vector &v);

    /// polygon[i] returns the vertex given by i mod number of vertices
    Vector operator[](int i) const;
    Vector& operator[](int i);

    /// Returns an outwards normal of the edge connecting vertex[i] to vertex[i+1]
    Vector GetOutwardsNormal(int i) const;

    /// Returns the directed area of the polygon
    double SignedArea() const;

    /// Returns (absolute) area of the polygon
    double AbsArea() const;

    /// Returns the centroid of the polygon (Assumes it is not self-intersecting)
    Vector Centroid() const;
};  

//// Free Functions on Polygons

// Returns the integral:
// Integral ||x - y||^2 dx
// With domain x inside the polygon P.
double EnergyIntegral(const Polygon &P, const Vector &y);
double OtherEnergyIntegral(const Polygon &P, const Vector &y);

Polygon CreateEncompassingSquare(const std::vector<Vector> &points);

Polygon CreateDiscretizedDisk(const Vector &x, const double r, const int num_sides);

// Shifts and scales the bounding box to be the [0,1]x[0,1] square. Then do the
// same transformation to all polygons in cells.
std::vector<Polygon> StandardizeCells(
        const std::vector<Polygon> &cells, 
        const Polygon &bounding_box);

/// Sutherland-Hodgman Polygon Clipping Algorithm
Polygon PolygonClip(Polygon subjectPolygon, Polygon clipPolygon);

// Diagram Algorithms
std::vector<Polygon> Voronoi(
        const std::vector<Vector> &voronoiSites,
        const Polygon &bounding_box);

std::vector<Polygon> PowerDiagram(
        const std::vector<Vector> &sites, 
        const std::vector<double> &weights,
        const Polygon &bounding_box);

std::vector<Polygon> FastPowerDiagram(
        const std::vector<Vector> &sites,
        const std::vector<double> &weights,
        const Polygon &bounding_box);

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
    );

    std::vector<Vector> LiftedKNN(Vector P, double weight, int k) const;

    std::vector<Polygon> MakeDiagram(const Polygon &bounding_box) const;
};

#endif