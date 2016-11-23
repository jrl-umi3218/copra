// This file is part of ModelPreviewController.

// ModelPreviewController is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ModelPreviewController is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ModelPreviewController.  If not, see
// <http://www.gnu.org/licenses/>.

#include "GUROBISolver.h"
#include <iostream>

namespace mpc {

/*
 * GUROBI
 */

GUROBISolver::GUROBISolver()
    : solver_(std::make_unique<Eigen::GurobiDense>())
{
}

int GUROBISolver::SI_fail() const
{
    return solver_->fail();
}

int GUROBISolver::SI_iter() const
{
    return solver_->iter();
}

void GUROBISolver::SI_inform() const
{
    std::cout << "Please add the inform() method" << std::endl;
}

const Eigen::VectorXd& GUROBISolver::SI_result() const
{
    return solver_->result();
}

void GUROBISolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool GUROBISolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& C,
    const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& Beq,
    const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& Bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    return solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq, XL, XU);
}

} // namespace pc