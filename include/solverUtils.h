// This file is part of copra.

// copra is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// copra is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with copra.  If not, see
// <http://www.gnu.org/licenses/>.

#pragma once

// stl
#include <memory>
#include <utility>

// copra
#include "QuadProgSolver.h"
#include "config.hh"

// optional mpc
#include "solverConfig.h"

namespace copra {

/**
 * Enum class that handles flag for selecting a qp solver.
 */
enum class SolverFlag {
    DEFAULT, /**< Default solver (QuadProgDense solver) */
#ifdef EIGEN_LSSOL_FOUND
    LSSOL, /**< Stanford LSSOL solver */
#endif
#ifdef EIGEN_GUROBI_FOUND
    GUROBIDense, /**< Gurobi quadratic dense solver */
#endif
#ifdef EIGEN_QLD_FOUND
    QLD, /**< Scilab QLD solver */
#endif
    QuadProgDense, /**< DenseMatrix version of QuadProg solver */
    // QuadProgSparse
};

/**
 * Helper function to get an unique pointer to a desired solver.
 * \param flag Flag of the solver.
 * \return An unique pointer to the desired solver.
 */
std::unique_ptr<SolverInterface> COPRA_DLLAPI solverFactory(SolverFlag flag);

/**
 * Helper function to get a desired solver.
 * This should only be used by python (unique_ptr are not bindable)
 * \param flag Flag of the solver.
 * \return The desired solver.
 */
SolverInterface* pythonSolverFactory(SolverFlag flag);

} // namespace pc