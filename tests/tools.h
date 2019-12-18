/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <Eigen/Core>
#include <utility>
#include <vector>

namespace tools {

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
