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

/**
 * Enum class that handles flag for selecting a qp solver.
 */
enum class SolverFlag
{
    DEFAULT, /**< Default solver (QuadProgDense solver) */
#ifdef LSSOL_SOLVER_FOUND
    LSSOL, /**< Standford LSSOL solver */
#endif
    QLD, /**< Scilab QLD solver */
    QuadProgDense, /**< DenseMatrix version of QuadProg solver */
    // QuadProgSparse
};

/**
 * Enum class that Allows selecting between two types of preview control.
 */
enum class PCFlag
{
    Full, /**< All the matrix is considered (from t=0 to t=tf) */
    Last /**< Only the last matrix is considered (at t=tf) */
};

/**
 * Helper function to get an unique pointer to a desired solver.
 * @param flag Flag of the solver. @see SolverFlag.
 * @return An unique pointer to the desired solver.
 */
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