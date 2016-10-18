#pragma once

#include "SolverInterface.h"
#include "QuadProgSolver.h"
#include "QLDSolver.h"
#ifdef LSSOL_SOLVER_FOUND
#include "LSSOLSolver.h"
#endif

#include <memory>
#include <utility>

namespace pc
{

enum class SolverFlag
{
    DEFAULT,
#ifdef LSSOL_SOLVER_FOUND
    LSSOL,
#endif
    QLD,
    QuadProgDense,
    // QuadProgSparse
};

enum class PCFlag
{
    Full,
    Last
};

std::unique_ptr<SolverInterface> solverFactory(SolverFlag flag)
{
    switch (flag)
    {
#ifdef LSSOL_SOLVER_FOUND
    case SolverFlag::LSSOL:
        return std::make_unique<LSSOLSolver>();
#endif
    case SolverFlag::QLD:
        return std::make_unique<QLDSolver>();
    // case SolverFlag::QuadProgSparse:
    //    return std::make_unique<QuadProgSparseSolver>();
    case SolverFlag::QuadProgDense:
    default:
        return std::make_unique<QuadProgDenseSolver>();
    }
}

} // namespace pc