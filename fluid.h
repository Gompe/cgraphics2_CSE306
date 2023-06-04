#pragma once

#ifndef FLUID_H
#define FLUID_H

#include "ot.h"

#include <vector>
#include <memory>

#include "gx_vector.h"
#include "gx_polygon.h"

class Fluid;

class FluidSolver : public ObjectiveLBFGS{
    Fluid *owner;
    int num_particles;

public:
    FluidSolver() = delete;
    FluidSolver(Fluid *owner);

    void InitializeVariablesArray();

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

    double ComputeFluidArea() const;
    double ComputeAirArea() const;

    double Functional() const;
    std::vector<double> GradientFunctional() const;

    void AdjustWeights();
    void GradientAscent(int n_iter, double lr);
};

class Fluid {
public:
    int num_particles;

    std::vector<Vector> positions;
    std::vector<Vector> velocities;

    Polygon bounding_box;

    double volume_fluid;
    double volume_air;

    // Physical constants
    const double k_min_diff_weights = 1E-3;

    const double k_volume_fluid = 0.3;
    const double k_volume_air = 0.7;
    const double k_elasticity = 0.5;
    const double k_epsilon = 0.01; // original 0.004
    const Vector k_g = Vector(0., -10., 0.); 

    std::unique_ptr<FluidSolver> solver;

    // To optimize
    std::vector<double> weights;
    double weight_air;

    std::vector<Polygon> cells;

    Fluid(const std::vector<Vector> &positions);
    Fluid(const std::vector<Vector> &positions, const Polygon &bounding_box);

    void InitializeParameters();

    void SaveFrame(const int frameid) const;

    void Run(const double dt, const int num_frames);
    
    // Physics
    void ComputeFluidCells();
    std::vector<Vector> ComputeForces() const;

private:
    void BounceParticlesBack();    
};

#endif