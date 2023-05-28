#include "ot.h"

namespace ot_detail {
    lbfgsfloatval_t Evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ) {
        return reinterpret_cast<ObjectiveLBFGS *>(instance)->Evaluate(x, g, n, step);
    }

    int Progress(
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
    ) {
        return reinterpret_cast<ObjectiveLBFGS *>(instance)->Progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }
}; // namespace ot_detail

ObjectiveLBFGS::~ObjectiveLBFGS() {
    FreeMemory();
}

void ObjectiveLBFGS::ReserveMemory(size_t num_variables) {
    if (capacity >= num_variables) {
        return;
    }

    FreeMemory();
    m_variables = lbfgs_malloc(num_variables);

    if (m_variables == nullptr) {
        fprintf(stderr, "LBFGS Malloc failed.\n");
        exit(1);
    }

    capacity = num_variables;
}

void ObjectiveLBFGS::FreeMemory() {
    if (m_variables != nullptr) {
        lbfgs_free(m_variables);
    }

    capacity = 0;
}

int ObjectiveLBFGS::CallLBFGS(int num_variables, lbfgsfloatval_t *fx_ptr) {
    return lbfgs(
        num_variables,
        m_variables,
        fx_ptr,
        ot_detail::Evaluate,
        ot_detail::Progress,
        this,
        nullptr
    );
}
