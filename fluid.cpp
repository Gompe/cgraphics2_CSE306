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

    ReserveMemory(num_particles + 1);

    InitializeVariablesArray();

    lbfgsfloatval_t fx;
    int ret = CallLBFGS(num_particles + 1, &fx);

    // printf("INFO: L-BFGS optimization terminated with status code = %d\n", ret);
    return ret;
}

void FluidSolver::InitializeVariablesArray() {
    for (int i = 0; i < num_particles; i++) {
        m_variables[i] = (lbfgsfloatval_t) owner->weights[i];
    }

    m_variables[num_particles] = (lbfgsfloatval_t) owner->weight_air;
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

    lbfgsfloatval_t fx = (lbfgsfloatval_t) Functional();
    std::vector<double> grad = GradientFunctional();

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

double FluidSolver::Functional() const {

    const double expected_particle_area = (owner->volume_fluid/(double) num_particles);

    // Terms coming from the fluid particles
    auto output_term = [this, expected_particle_area](auto i){
        return EnergyIntegral(this->owner->cells[i], this->owner->positions[i]) +
                  this->owner->weights[i] * (expected_particle_area - this->owner->cells[i].AbsArea());
    };

    double output = 0.;
    #pragma omp parallel for reduction(+:output)
    for (size_t i = 0; i < num_particles; i++) {
        const auto term = output_term(i);
        output += term;
    }

    // for (int i = 0; i < num_particles; i++) {
    //     output += EnergyIntegral(owner->cells[i], owner->positions[i]) +
    //               owner->weights[i] * (expected_particle_area - owner->cells[i].AbsArea());
    // }

    // Terms coming from the air particles
    output += owner->weight_air * (owner->volume_air - ComputeAirArea());

    return output;
}

std::vector<double> FluidSolver::GradientFunctional() const {
    std::vector<double> gradient(num_particles + 1);

    for (int i = 0; i < num_particles; i++) {
        gradient[i] = (owner->volume_fluid/(double) num_particles) - owner->cells[i].AbsArea();
    }

    // Air weight
    gradient[num_particles] = owner->volume_air - ComputeAirArea();
    return gradient;
}

void FluidSolver::AdjustWeights() {
    const double weight_air = m_variables[num_particles];
    const double weight_lb = weight_air + (owner->volume_fluid/(double) num_particles)/M_PI;

    auto adjust_weight = [&](lbfgsfloatval_t &w){
        w = std::max(w, (lbfgsfloatval_t) weight_lb) - (lbfgsfloatval_t) weight_air;
    };
    
    std::for_each(m_variables, m_variables + num_particles, adjust_weight);
    m_variables[num_particles] = 0.;
}

void FluidSolver::GradientAscent(int n_iter, double lr) {
    ReserveMemory(num_particles + 1);

    for (int k = 0; k < num_particles; k++) {
        m_variables[k] = owner->weights[k];
    }
    m_variables[num_particles] = owner->weight_air;

    // double fx = Functional();
    // std::cout << "GA enter: fx = " << fx << std::endl;

    for (int iter = 0; iter < n_iter; iter++) {
        std::vector<double> grad = GradientFunctional();
        // Not updating w_air
        for (int k = 0; k < num_particles; k++) {
            m_variables[k] += (lbfgsfloatval_t) lr * grad[k];
        }
        // Make sure that weight_air is not larger than any fluid weight
        // AdjustWeights();

        UpdateWeights();
        UpdateCells();
    }

    // fx = Functional();
    // std::cout << "GA leave: fx = " << fx << std::endl;
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

    const double box_area = bounding_box.AbsArea();
    volume_fluid = k_volume_fluid * box_area;
    volume_air = k_volume_air * box_area;

    const double estimated_cell_area = (volume_fluid/(double) num_particles) / M_PI;

    weights.resize(num_particles);
    std::fill(weights.begin(), weights.end(), estimated_cell_area);

    weight_air = 0;

    cells.resize(num_particles);
    solver.reset(new FluidSolver(this)); 
}

void Fluid::SaveFrame(const int frameid) const {
    std::vector<Polygon> particles;

    particles.insert(particles.begin(), cells.begin(), cells.end());
    particles = StandardizeCells(particles, bounding_box);

    save_frame(particles, "./frames/frame", frameid);
    std::cout << "Saved frame " << frameid << std::endl;
}

void Fluid::Run(const double dt, const int num_frames) {

    auto time_step_update = [dt](auto &X, auto &X_prime) {
        auto update_one = [dt](auto x, auto x_prime) {
            return x + dt * x_prime;
        };

        std::transform(X.begin(), X.end(), X_prime.begin(), X.begin(),
                       update_one);
    };

    for (int frameid = 0; frameid < num_frames; frameid++) {
        ComputeFluidCells();

        // Show the particles and the cells
        SaveFrame(frameid); 

        std::vector<Vector> forces = ComputeForces();
        
        const double mass = 200 * (volume_fluid / (double) num_particles);
        std::vector<Vector> accelerations(num_particles);

        std::transform(forces.begin(), forces.end(), accelerations.begin(), 
                       [mass](auto f){ return f / mass; });

        // Update velocites
        time_step_update(velocities, accelerations);

        //Update Positions
        time_step_update(positions, velocities);

        // Handle particles out of the domain
        BounceParticlesBack();
    }
}

void Fluid::ComputeFluidCells() {
    // Run Optimal Transport Algorithm
    solver->Run();

    // Fine-tune with gradient descent
    solver->GradientAscent(100, 0.05);
}

std::vector<Vector> Fluid::ComputeForces() const {
    const double mass = 200; // * (volume_fluid / (double) num_particles);
    std::vector<Vector> forces;

    for (int i = 0; i < num_particles; i++) {
        const Vector centroid = [&](){
            if (this->cells[i].size() == 0) {
                return this->positions[i];
            }
            else {
                return this->cells[i].Centroid();
            }
        }();

        forces.push_back(mass * k_g + (centroid - positions[i])/(k_epsilon * k_epsilon));
    }

    return forces;
}

void Fluid::BounceParticlesBack() {
    double min_x = polygon_detail::MinX(bounding_box.vertices);
    double min_y = polygon_detail::MinY(bounding_box.vertices);
    double max_x = polygon_detail::MaxX(bounding_box.vertices);
    double max_y = polygon_detail::MaxY(bounding_box.vertices);

    for (int i = 0; i < num_particles; i++) {
        // Left Border
        if (positions[i][0] < min_x) {
            positions[i][0] = min_x + k_elasticity * (min_x - positions[i][0]);
            velocities[i][0] *= -k_elasticity;
        }
        // Right Border
        if (positions[i][0] > max_x) {
            positions[i][0] = max_x + k_elasticity * (max_x - positions[i][0]);
            velocities[i][0] *= -k_elasticity;
        }

        // Bottom Border
        if (positions[i][1] < min_y) {
            positions[i][1] = min_y + k_elasticity * (min_y - positions[i][1]);
            velocities[i][1] *= -k_elasticity;
        }

        // Top Border
        if (positions[i][1] > max_y) {
            positions[i][1] = max_y + k_elasticity * (max_y - positions[i][1]);
            velocities[i][1] *= -k_elasticity;
        }
    }
}