#pragma once

#ifndef LAGUERRE_H
#define LAGUERRE_H

#include "ot.h"

#include <memory>

#include "gx_vector.h"
#include "gx_polygon.h"



class LaguerreDiagram;

class LaguerreSolver : public ObjectiveLBFGS {
    LaguerreDiagram *owner;
    int num_sites;

public:
    LaguerreSolver() = delete;
    LaguerreSolver(LaguerreDiagram *owner);
    int Run();

    lbfgsfloatval_t Evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    );

    int Progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    );

    void UpdateWeights();
    void UpdateCells();

    double Functional() const;
    std::vector<double> GradientFunctional() const;
};

class LaguerreDiagram {
    std::unique_ptr<LaguerreSolver> solver;
public:
    std::vector<Vector> sites;
    std::vector<Polygon> cells;
    std::vector<double> lambdas;
    std::vector<double> weights;
    Polygon bounding_box;

    LaguerreDiagram(
        const std::vector<Vector> &sites,
        const std::vector<Polygon> &cells,
        const std::vector<double> &lambdas,
        const std::vector<double> &weights
    );

    LaguerreDiagram(
        const std::vector<Vector> &sites,
        const std::vector<Polygon> &cells,
        const std::vector<double> &lambdas,
        const std::vector<double> &weights,
        const Polygon &bounding_box
    );

    int Optimize();
};

#endif