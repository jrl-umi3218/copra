// This file is part of mpc.

// mpc is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with mpc.  If not, see
// <http://www.gnu.org/licenses/>.
#pragma once

// stl
#include <memory>
#include <utility>

// mpc
#include "config.hh"
#include "QuadProgSolver.h"

// optional mpc
#ifdef QLD_SOLVER_FOUND
#include "QLDSolver.h"
#endif
#include "SolverInterface.h"
#ifdef LSSOL_SOLVER_FOUND
#include "LSSOLSolver.h"
#endif
#ifdef GUROBI_SOLVER_FOUND
#include "GUROBISolver.h"
#endif


namespace mpc {

/**
 * Enum class that handles flag for selecting a qp solver.
 */
enum class SolverFlag {
    DEFAULT, /**< Default solver (QuadProgDense solver) */
#ifdef LSSOL_SOLVER_FOUND
    LSSOL, /**< Standford LSSOL solver */
#endif
#ifdef GUROBI_SOLVER_FOUND
    GUROBIDense, /**< Gurobi quadratic dense solver */
#endif
#ifdef QLD_SOLVER_FOUND
    QLD, /**< Scilab QLD solver */
#endif
    QuadProgDense, /**< DenseMatrix version of QuadProg solver */
    // QuadProgSparse
};

/**
 * Helper function to get an unique pointer to a desired solver.
 * @param flag Flag of the solver.
 * @return An unique pointer to the desired solver.
 */
std::unique_ptr<SolverInterface> MPC_DLLAPI solverFactory(SolverFlag flag)
{
    switch (flag) {
#ifdef LSSOL_SOLVER_FOUND
    case SolverFlag::LSSOL:
        return std::make_unique<LSSOLSolver>();
#endif
#ifdef GUROBI_SOLVER_FOUND
    case SolverFlag::GUROBIDense:
        return std::make_unique<GUROBISolver>();
#endif
#ifdef QLD_SOLVER_FOUND
    case SolverFlag::QLD:
        return std::make_unique<QLDSolver>();
#endif
    // case SolverFlag::QuadProgSparse:
    //    return std::make_unique<QuadProgSparseSolver>();
    case SolverFlag::QuadProgDense:
    default:
        return std::make_unique<QuadProgDenseSolver>();
    }
}

/**
 * Helper function to get a desired solver.
 * This should only be used by python (unique_ptr are not yet bindable)
 * @param flag Flag of the solver.
 * @return The desired solver.
 */
SolverInterface* pythonSolverFactory(SolverFlag flag)
{
    switch (flag) {
#ifdef LSSOL_SOLVER_FOUND
    case SolverFlag::LSSOL:
        return new LSSOLSolver;
#endif
#ifdef GUROBI_SOLVER_FOUND
    case SolverFlag::GUROBIDense:
        return new GUROBISolver;
#endif
#ifdef QLD_SOLVER_FOUND
    case SolverFlag::QLD:
        return new QLDSolver;
#endif
    case SolverFlag::QuadProgDense:
    default:
        return new QuadProgDenseSolver;
    }
}

} // namespace pc