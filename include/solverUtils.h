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

#include "api.h"

#include "QuadProgSolver.h"
#ifdef EIGEN_QLD_FOUND
#include "QLDSolver.h"
#endif
#ifdef EIGEN_LSSOL_FOUND
#include "LSSOLSolver.h"
#endif
#ifdef EIGEN_GUROBI_FOUND
#include "GUROBISolver.h"
#endif
#include <memory>
#include <utility>

namespace copra {

/**
 * Enum class that handles flag for selecting a qp solver.
 */
enum class COPRA_DLLAPI SolverFlag {
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
COPRA_DLLAPI std::unique_ptr<SolverInterface> solverFactory(SolverFlag flag);

/**
 * Helper function to get a desired solver.
 * This should only be used by python (unique_ptr are not bindable)
 * \param flag Flag of the solver.
 * \return The desired solver.
 */
COPRA_DLLAPI SolverInterface* pythonSolverFactory(SolverFlag flag);

} // namespace pc