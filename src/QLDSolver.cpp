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

#include "QLDSolver.h"
#include <iostream>

namespace copra {

/*
 * QLD
 */

QLDSolver::QLDSolver()
    : solver_(new Eigen::QLD())
    , eps_(1e-8)
{
}

int QLDSolver::SI_fail() const
{
    return solver_->fail();
}

void QLDSolver::SI_inform() const
{
    switch (solver_->fail()) {
    case 0:
        std::cout << "The optimality conditions are satisfied." << std::endl;
        break;
    case 1:
        std::cout << "The algorithm has been stopped after too many "
                     "MAXIT iterations (40*(N+M))."
                  << std::endl;
        break;
    case 2:
        std::cout << "Termination accuracy insufficient to satisfy "
                     "convergence criterion."
                  << std::endl;
        break;
    case 3:
        std::cout << "Internal inconsistency of QL, division by zero." << std::endl;
        break;
    case 4:
        std::cout << "Numerical instability prevents successful termination. "
                     "Use tolerance specified in WAR(1) for a restart."
                  << std::endl;
        break;
    case 5:
        std::cout << "Length of a working array is too short." << std::endl;
        break;
    default:
        std::cout << "Constraints are inconsistent and IFAIL=100+ICON, "
                     "where ICON denotes a constraint causing the conflict."
                  << std::endl;
        break;
    }
}

void QLDSolver::SI_printLevel(int pl)
{
    solver_->verbose(pl != 0);
}

void QLDSolver::SI_feasibilityTolerance(double tol)
{
    eps_ = tol;
}

const Eigen::VectorXd&
QLDSolver::SI_result() const
{
    return solver_->result();
}

void QLDSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool QLDSolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
    const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
    const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    return solver_->solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU, eps_);
}

} // namespace pc