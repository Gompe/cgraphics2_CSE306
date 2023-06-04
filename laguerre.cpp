#include "laguerre.h"

LaguerreSolver::LaguerreSolver(LaguerreDiagram *owner) : 
    owner(owner) 
{
    num_sites = owner->sites.size();
}

int LaguerreSolver::Run() {
    lbfgsfloatval_t fx;

    ReserveMemory(num_sites);

    // Initialize variables array
    for (int i = 0; i < num_sites; i++) {
        m_variables[i] = (lbfgsfloatval_t) owner->weights[i];
    }

    int ret = CallLBFGS(num_sites, &fx);

    // printf("INFO: L-BFGS optimization terminated with status code = %d\n", ret);

    return ret;
}

lbfgsfloatval_t LaguerreSolver::Evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step) {
    
    assert(x == m_variables);
    assert(n == num_sites);

    UpdateWeights();
    UpdateCells();

    lbfgsfloatval_t fx = 0.0;

    fx = (lbfgsfloatval_t) Functional();
    std::vector<double> grad = GradientFunctional();

    // Swap signs because we are maximizing 
    fx = -fx;
    for (int i = 0; i < num_sites; i++) {
        g[i] = (lbfgsfloatval_t) -grad[i];
    }

    return fx;
}

void LaguerreSolver::UpdateWeights() {
    for (int i=0; i < num_sites; i++) {
        owner->weights[i] = (double) m_variables[i];
    }
}

void LaguerreSolver::UpdateCells() {
    owner->cells = FastPowerDiagram(owner->sites, owner->weights, owner->bounding_box);
}

double LaguerreSolver::Functional() const {
    double output = 0.;

    double totalArea = std::accumulate(owner->cells.begin(), owner->cells.end(), 0.,
        [](double abs, const Polygon& cell) {return abs + cell.AbsArea(); });

    for (int i = 0; i < num_sites; i++) {
        const Polygon &cell = owner->cells[i]; // this cell
        double integral = EnergyIntegral(cell, owner->sites[i]);
        integral -= cell.AbsArea() * owner->weights[i];

        // Adding to output
        output += integral;
        output += owner->lambdas[i] * owner->weights[i];
    }

    return output;
}

std::vector<double> LaguerreSolver::GradientFunctional() const {
    std::vector<double> output(num_sites);

    double totalArea = std::accumulate(owner->cells.begin(), owner->cells.end(), 0.,
        [](double abs, const Polygon& cell) {return abs + cell.AbsArea(); });

    for (int i = 0; i < num_sites; i++) {
        output[i] = owner->lambdas[i] - owner->cells[i].AbsArea();
    }

    return output;
}

int LaguerreSolver::Progress(
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
/////////////// Laguerre Diagram Implementation ///////////////////////////////

LaguerreDiagram::LaguerreDiagram(
        const std::vector<Vector> &sites,
        const std::vector<Polygon> &cells,
        const std::vector<double> &lambdas,
        const std::vector<double> &weights
    ) :
    sites(sites),
    cells(cells),
    lambdas(lambdas),
    weights(weights)
{
    this->bounding_box = CreateEncompassingSquare(sites);
    solver.reset(new LaguerreSolver(this));
}

LaguerreDiagram::LaguerreDiagram(
        const std::vector<Vector> &sites,
        const std::vector<Polygon> &cells,
        const std::vector<double> &lambdas,
        const std::vector<double> &weights,
        const Polygon &bounding_box
    ) :
    sites(sites),
    cells(cells),
    lambdas(lambdas),
    weights(weights),
    bounding_box(bounding_box)
{
    solver.reset(new LaguerreSolver(this));
}

int LaguerreDiagram::Optimize() {
    return solver->Run();
}