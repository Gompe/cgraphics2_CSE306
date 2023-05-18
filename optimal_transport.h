#include "lib_lbfgs/lbfgs.h"
#include "gx_polygon.h"

#include <iostream>

#ifndef LBFGS_VERBOSE
    #define LBFGS_VERBOSE 0
#endif

namespace optimal_transport {

static double FunctionalG(
    const std::vector<Vector> &sites,
    const std::vector<Polygon> &cells,
    const std::vector<double> &lambdas,
    const std::vector<double> &weights
    )
{
    double output = 0.;
    const int N = sites.size();

    double totalArea = std::accumulate(cells.begin(), cells.end(), 0.,
        [](double abs, const Polygon& cell) {return abs + cell.absArea(); });

    for (int i=0; i<N; i++) {
        // Pages 102, 103 of lecture notes
        const Polygon &cell = cells[i]; // this cell

        // Computing integral
        double integral = 0.;
        for (int k=0; k < cell.size(); k++) {
            double leftFactor = (cell[k-1][0]*cell[k][1] - cell[k][0]*cell[k-1][1]);
            double rightFactor = cell[k-1].norm2() + cell[k].norm2() + dot(cell[k-1], cell[k]);

            rightFactor += -4 * dot(sites[i], cell[k-1] + cell[k]) + 6 * sites[i].norm2();

            integral += leftFactor * rightFactor;
        } 
        integral /= 12.;
        integral = std::abs(integral);
        integral -= cell.absArea() * weights[i];

        // Scale integral by total Area
        integral /= totalArea;

        // Adding to output
        output += integral;
        output += lambdas[i] * weights[i];
    }

    return output;
}

static std::vector<double> GradFunctionalG(
    const std::vector<Vector> &sites,
    const std::vector<Polygon> &cells,
    const std::vector<double> &lambdas,
    const std::vector<double> &weights
    )
{
    const int N = sites.size();
    std::vector<double> output(N);

    double totalArea = std::accumulate(cells.begin(), cells.end(), 0.,
        [](double abs, const Polygon& cell) {return abs + cell.absArea(); });

    for (int i=0; i<N; i++) {
        output[i] = lambdas[i] - cells[i].absArea()/totalArea;
    }

    return output;
}

class ObjectiveFunction
{
    lbfgsfloatval_t *m_weightArray = nullptr;

    const int N;
    const std::vector<Vector> sites;
    std::vector<Polygon> &cells;
    const std::vector<double> lambdas;
    std::vector<double> &weights;

    void m_lbfgs_allocate(){
        if (m_weightArray == nullptr)
            m_weightArray = lbfgs_malloc(N);
        if (m_weightArray == nullptr) {
            fprintf(stderr, "ERROR: Failed lbfgs_malloc.\n");
            exit(1);
        }
    }

    void m_lbfgs_deallocate(){
        if (m_weightArray != nullptr) 
                    lbfgs_free(m_weightArray);
    }

public:
    ObjectiveFunction(
        const std::vector<Vector> &sites,
        std::vector<Polygon> &cells,
        const std::vector<double> &lambdas,
        std::vector<double> &weights
    ) : sites(sites), cells(cells), lambdas(lambdas), weights(weights), N(sites.size())
    {}

    ~ObjectiveFunction()
    {
        m_lbfgs_deallocate();
    }

    int run()
    {
        lbfgsfloatval_t fx;

        m_lbfgs_allocate();

        // Initialize array
        for (int i=0; i<N; i++) {
            m_weightArray[i] = (lbfgsfloatval_t) weights[i];
        }

        int ret = lbfgs(N, m_weightArray, &fx, _evaluate, _progress, this, nullptr);
        // m_UpdateParams();

        #if LBFGS_VERBOSE
            printf("INFO: L-BFGS optimization terminated with status code = %d\n", ret);
        #endif

        return ret;
    }
private:
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<ObjectiveFunction*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        assert(x == m_weightArray);
        assert(n == N);

        m_UpdateParams();

        lbfgsfloatval_t fx = 0.0;

        // Don't forget to swap the sign for both G and GradG
        fx = (lbfgsfloatval_t) -FunctionalG(sites, cells, lambdas, weights);
        std::vector<double> grad = GradFunctionalG(sites, cells, lambdas, weights);

        for (int i=0; i<N; i++) {
            g[i] = (lbfgsfloatval_t) -grad[i];
        }

        return fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<ObjectiveFunction*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {   
        #if LBFGS_VERBOSE
            printf("Iteration %d:\n", k);
            printf("  fx = %f\n", fx);
            printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
            printf("\n");
        #endif

        return 0;
    }

    void m_UpdateWeights()
    {
        for (int i=0; i<N; i++) {
            weights[i] = (double) m_weightArray[i];
        }

        double minWeight = *std::min_element(weights.begin(), weights.end());
        
        std::for_each(weights.begin(), weights.end(), 
            [minWeight](double &weight){ weight -= minWeight; });
    }

    void m_UpdateCells()
    {
        cells = powerDiagram(sites, weights);
    }

    void m_UpdateParams()
    {
        m_UpdateWeights();
        m_UpdateCells();
    }

};  

}; // namespace optimal_transport

extern int optimalTransport(
    const std::vector<Vector> &sites,
    std::vector<Polygon> &cells,
    const std::vector<double> &lambdas,
    std::vector<double> &weights
    )
{
    auto obj = optimal_transport::ObjectiveFunction(
        sites,
        cells,
        lambdas,
        weights
    );

    return obj.run();
}