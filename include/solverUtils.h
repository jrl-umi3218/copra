/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

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
#ifdef EIGEN_OSQP_FOUND
#include "OSQPSolver.h"
#endif
#include <memory>
#include <utility>

namespace copra {

#ifdef __GNUC__
#pragma GCC diagnostic push
// Work around GCC (< 6) bug see: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=43407
#pragma GCC diagnostic ignored "-Wattributes"
#endif

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
#ifdef EIGEN_OSQP_FOUND
    OSQP,
#endif
    QuadProgDense, /**< DenseMatrix version of QuadProg solver */
    // QuadProgSparse
};

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

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
