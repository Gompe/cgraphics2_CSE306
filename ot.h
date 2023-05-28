#pragma once

#ifndef OT_H
#define OT_H

#include "lib_lbfgs/lbfgs.h"

#include <stdlib.h>

#include <cstddef>
#include <cstdio>

namespace ot_detail {
    lbfgsfloatval_t Evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    );
    
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
    );
}; // namespace ot_detail

/// Abstract Class For Objective Functions
class ObjectiveLBFGS {
protected:
    // Array containing variables to be optimized by the LBFGS engine
    lbfgsfloatval_t *m_variables = nullptr;
    size_t capacity = 0;

    void ReserveMemory(size_t num_variables);
    void FreeMemory();

public:
    ObjectiveLBFGS() = default;
    ~ObjectiveLBFGS();

    // Runs the optimization procedure
    virtual int Run() = 0;
    int CallLBFGS(int num_variables, lbfgsfloatval_t *fx);

    virtual lbfgsfloatval_t Evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
    ) = 0;

    virtual int Progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
    ) = 0;
};

#endif