/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "solverUtils.h"
#include <Eigen/Core>
#include <utility>
#include <vector>

namespace tools {

static const std::vector<std::pair<std::string, copra::SolverFlag>> Solvers = {
    { "Default (QuadProgDense)", copra::SolverFlag::DEFAULT },
#ifdef EIGEN_LSSOL_FOUND
    { "LSSOL", copra::SolverFlag::LSSOL },
#endif
#ifdef EIGEN_QLD_FOUND
    { "QLD", copra::SolverFlag::QLD },
#endif
#ifdef EIGEN_GUROBI_FOUND
    { "GUROBIDense", copra::SolverFlag::GUROBIDense },
#endif
#ifdef EIGEN_OSQP_FOUND
    { "OSQP", copra::SolverFlag::OSQP },
#endif
};

using solver_timers_t = std::vector<std::pair<std::string, double>>;

struct SolverTimers {
    solver_timers_t st; // solve time
    solver_timers_t bt; // build time
    solver_timers_t ct; // coombine time
};

Eigen::MatrixXd spanMatrix(const Eigen::MatrixXd& m, int size, int addCols = 0);

Eigen::VectorXd spanVector(const Eigen::VectorXd& v, int size);

std::string getSortedTimers(SolverTimers& solT);

} // namespace
