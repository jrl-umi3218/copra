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

#include "LSSOLSolver.h"

namespace mpc {

/*
 * LSSOL
 */

LSSOLSolver::LSSOLSolver()
    : solver_(std::make_unique<Eigen::StdLSSOL>())
{
}

int LSSOLSolver::SI_fail() const
{
    return solver_->fail();
}

int LSSOLSolver::SI_iter() const
{
    return solver_->iter();
}

void LSSOLSolver::SI_inform() const
{
    solver_->inform();
}

void LSSOLSolver::SI_printLevel(int pl)
{
    solver_->printLevel(pl);
}

void LSSOLSolver::SI_feasibilityTolerance(double tol)
{
    solver_->feasibilityTol(tol);
}

bool LSSOLSolver::SI_warmStart() const
{
    return solver_->warm();
}

void LSSOLSolver::SI_warmStart(bool w)
{
    solver_->warm(w);
}

const Eigen::VectorXd&
LSSOLSolver::SI_result() const
{
    return solver_->result();
}

void LSSOLSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool LSSOLSolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
    const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
    const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    return solver_->solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU);
}

} // namespace pc