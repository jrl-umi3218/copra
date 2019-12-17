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

#include "OSQPSolver.h"
#include <iostream>

namespace copra
{

/*
 * OSQP
 */

OSQPSolver::OSQPSolver()
{
}

int OSQPSolver::SI_fail() const
{
    return static_cast<int>(solver_.status());
}

void OSQPSolver::SI_inform() const
{
    solver_.inform(std::cout);
}

void OSQPSolver::SI_printLevel(int pl)
{
    solver_.verbose(pl != 0);
}

void OSQPSolver::SI_feasibilityTolerance(double tol)
{
    solver_.absConvergenceTol(tol);
}

const Eigen::VectorXd& OSQPSolver::SI_result() const
{
    return solver_.result();
}

void OSQPSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_.problem(nrVar, nrEq + nrInEq);
}

bool OSQPSolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
    const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
    const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    Eigen::MatrixXd A(Aeq.rows() + Aineq.rows(), Aeq.cols());
    A.topRows(Aeq.rows()) = Aeq;
    A.bottomRows(Aineq.rows()) = Aineq;
    Eigen::VectorXd bl(Aeq.rows() + Aineq.rows());
    Eigen::VectorXd bu(Aeq.rows() + Aineq.rows());
    bl.head(Aeq.rows()) = beq;
    bl.tail(Aineq.rows()).setConstant(-std::numeric_limits<double>::max());
    bu.head(Aeq.rows()) = beq;
    bu.tail(Aineq.rows()) = bineq;

    return solver_.solve(Q, c, A, bl, bu, XL, XU);
}

} // namespace copra
