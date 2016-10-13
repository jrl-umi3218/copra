#include "SolverInterface.h"

namespace pc {

enum class SolverFlag {
    DEFAULT,
#ifdef LSSOL_SOLVER_FOUND
    LSSOL,
#endif
    QLD,
    QuadProgDense,
    // QuadProgSparse
};

enum class PCFlag { Full,
    Last };

static const std::unique_ptr<SolverInterface> solverFactory(SolverFlag)
{
    switch (flag) {
    case SolverFlag::DEFAULT:
        return std::make_unique<QuadProgDenseSolver>();
#ifdef LSSOL_SOLVER_FOUND
    case SolverFlag::LSSOL:
        return std::make_unique<LSSOLSolver>();
#endif
    case SolverFlag::QLD:
        return std::make_unique<QLDSolver>();
    case SolverFlag::QuadProgDense:
        return std::make_unique<QuadProgDenseSolver>();
        // case SolverFlag::QuadProgSparse:
        //    return std::make_unique<QuadProgSparseSolver>();
    }
}
}