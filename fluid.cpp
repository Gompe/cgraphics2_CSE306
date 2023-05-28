#include "fluid.h"
#include "svg_handler.h"

#include <algorithm>

#define GAMMA 1

FluidSolver::FluidSolver(Fluid *owner) :
    owner(owner)
{
    num_particles = owner->num_particles;
}

int FluidSolver::Run() {
    lbfgsfloatval_t fx;

    ReserveMemory(num_particles + 1);

    // Initialize variables array
    for (int i = 0; i < num_particles; i++) {
        m_variables[i] = (lbfgsfloatval_t) owner->weights[i]/owner->bounding_box.AbsArea();
    }

    m_variables[num_particles] = owner->weight_air/owner->bounding_box.AbsArea();

    int ret = CallLBFGS(num_particles + 1, &fx);

    printf("INFO: L-BFGS optimization terminated with status code = %d\n", ret);
    return ret;
}

lbfgsfloatval_t FluidSolver::Evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step) {
    
    assert(x == m_variables);
    assert(n == num_particles + 1);

    UpdateWeights();
    UpdateCells();

    lbfgsfloatval_t fx = 0.0;

    fx = (lbfgsfloatval_t) Functional(x);
    std::vector<double> grad = GradientFunctional(x);

    // Swap signs because we are maximizing 
    fx = -fx;
    for (int i = 0; i < num_particles + 1; i++) {
        g[i] = (lbfgsfloatval_t) -grad[i];
    }

    return fx;
}

void FluidSolver::UpdateWeights() {
    for (int i=0; i < num_particles; i++) {
        owner->weights[i] = (double) m_variables[i];
    }
    owner->weight_air = (double) m_variables[num_particles];
}

static std::vector<Polygon> AirPowerDiagram(
        const std::vector<Vector> &sites,
        const std::vector<double> &weights,
        const double weight_air,
        const Polygon &bounding_box) {
    
    std::vector<Polygon> laguerre_cells = FastPowerDiagram(sites, weights, bounding_box);

    const int N = sites.size();
    for (int i = 0; i < N; i++) {
        if (weights[i] <= weight_air) {
            laguerre_cells[i] = Polygon();
        }
        else {
            Polygon disk = CreateDiscretizedDisk(sites[i], sqrt(weights[i] - weight_air), polygon_detail::num_sides_circle);
            laguerre_cells[i] = PolygonClip(laguerre_cells[i], disk);
        }
    }

    return laguerre_cells;
}

void FluidSolver::UpdateCells() {
    owner->cells = AirPowerDiagram(
        owner->positions, 
        owner->weights, 
        owner->weight_air,
        owner->bounding_box);
}

double FluidSolver::ComputeFluidArea() const {
    double fluid_area = 0.;
    for (const auto &cell : owner->cells) {
        fluid_area += cell.AbsArea();
    }

    return fluid_area;
}

double FluidSolver::ComputeAirArea() const {
    double box_area = owner->bounding_box.AbsArea();
    return box_area - ComputeFluidArea();
}

double FluidSolver::Functional(const lbfgsfloatval_t *W) const {
    double output = 0.;

    // Terms coming from the fluid particles
    for (int i = 0; i < num_particles; i++) {
        // Computing integral
        double integral = EnergyIntegral(owner->cells[i], owner->positions[i]);
        integral -= owner->cells[i].AbsArea() * owner->weights[i];

        // Adding to output
        output += integral;
        output += (owner->volume_fluid/num_particles) * owner->weights[i];
    }

    // Terms coming from the air particles
    output += owner->volume_air * owner->weight_air;
    output -= owner->weight_air * ComputeAirArea();

    return output;
}

std::vector<double> FluidSolver::GradientFunctional(const lbfgsfloatval_t *W) const {
    std::vector<double> gradient(num_particles + 1);

    for (int i = 0; i < num_particles; i++) {
        gradient[i] = (owner->volume_fluid/num_particles) - owner->cells[i].AbsArea();
    }

    // Air weight
    gradient[num_particles] = owner->volume_air - ComputeAirArea();
    return gradient;
}

void FluidSolver::GradientAscent(int n_iter, double lr) {
    ReserveMemory(num_particles + 1);

    for (int k = 0; k < num_particles; k++) {
        m_variables[k] = owner->weights[k]/owner->bounding_box.AbsArea();
    }
    m_variables[num_particles] = owner->weight_air/owner->bounding_box.AbsArea();

    UpdateWeights();
    UpdateCells();

    double fx = Functional(m_variables);

    for (int iter = 0; iter < n_iter; iter++) {
        std::vector<double> grad = GradientFunctional(m_variables);
        for (int k = 0; k < num_particles + 1; k++) {
            m_variables[k] += (lbfgsfloatval_t) lr * grad[k];
        }
        UpdateWeights();
        UpdateCells();
    }

    fx = Functional(m_variables);
}

int FluidSolver::Progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls) 
{
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//// Fluid Implementation

Fluid::Fluid(const std::vector<Vector> &positions) :
        positions(positions) {
    this->bounding_box = CreateEncompassingSquare(positions);
    InitializeParameters();
}

Fluid::Fluid(const std::vector<Vector> &positions, const Polygon &bounding_box) :
        positions(positions),
        bounding_box(bounding_box) {
    InitializeParameters();
}

void Fluid::InitializeParameters() {
    num_particles = positions.size();

    velocities.resize(num_particles);

    weights.resize(num_particles);
    std::fill(weights.begin(), weights.end(), 1);

    weight_air = 0;

    double box_area = bounding_box.AbsArea();
    volume_fluid = k_volume_fluid * box_area;
    volume_air = k_volume_air * box_area;

    cells.resize(num_particles);
    solver.reset(new FluidSolver(this)); 
}

void Fluid::Run(double dt, int num_frames) {

    // std::vector<Polygon> particles;
    // for (int i = 0; i < num_particles; i++) {
    //     particles.push_back(CreateDiscretizedDisk(positions[i], 0.05, 50));
    // }
    // particles = StandardizeCells(particles, bounding_box);
    // save_frame(particles, "./frames/frame", 0);

    for (int frameid = 1; frameid < num_frames; frameid++) {
        OneStep(dt);

        std::vector<Polygon> particles;

        for (int i = 0; i < num_particles; i++) {
            particles.push_back(CreateDiscretizedDisk(positions[i], 0.01, 50));
        }
        particles.insert(particles.begin(), cells.begin(), cells.end());
        particles = StandardizeCells(particles, bounding_box);

        save_frame(particles, "./frames/frame", frameid);
    }
}

void Fluid::OneStep(double dt) {
    // weights, weight_air = OT

    weight_air = 0;
    std::fill(weights.begin(), weights.end(), 1);
    
    solver->Run();
    // solver->GradientAscent(1000, 1E-8);

    std::vector<double> DEBUG_ACCELERATIONS;

    for (int i = 0; i < num_particles; i++) {
        // Mass of the particle
        double mass = 200 * (volume_fluid / (double) num_particles);

        // Force
        Vector centroid;

        if (cells[i].size() == 0) {
            centroid = Vector();
        }
        else {
            centroid = cells[i].Centroid();
        }

        Vector F = 1/(k_epsilon * k_epsilon) * (centroid - positions[i]);
        F += mass * k_g;

        // Get acceleration (Newton's Law)
        Vector a = F/mass;
        velocities[i] = velocities[i] + dt * a;

        // Update Position
        positions[i] = positions[i] + dt * velocities[i];
    }

    this->BounceParticlesBack();
}

void Fluid::BounceParticlesBack() {
    double min_x = polygon_detail::MinX(bounding_box.vertices);
    double min_y = polygon_detail::MinY(bounding_box.vertices);
    double max_x = polygon_detail::MaxX(bounding_box.vertices);
    double max_y = polygon_detail::MaxY(bounding_box.vertices);

    for (int i = 0; i < num_particles; i++) {
        // Left Border
        if (positions[i][0] < min_x) {
            positions[i][0] = min_x + (min_x - positions[i][0]);
            velocities[i][0] *= -k_elasticity;
        }
        // Right Border
        if (positions[i][0] > max_x) {
            positions[i][0] = max_x + (max_x - positions[i][0]);
            velocities[i][0] *= -k_elasticity;
        }

        // Bottom Border
        if (positions[i][1] < min_y) {
            positions[i][1] = min_y + (min_y - positions[i][1]);
            velocities[i][1] *= -k_elasticity;
        }

        // Top Border
        if (positions[i][1] > max_y) {
            positions[i][1] = max_y + (max_y - positions[i][1]);
            velocities[i][1] *= -k_elasticity;
        }
    }
}