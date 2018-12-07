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

#include "LSSOLSolver.h"

namespace copra {

/*
 * LSSOL
 */

LSSOLSolver::LSSOLSolver()
    : solver_(new Eigen::LSSOL_QP())
{
}

int LSSOLSolver::SI_fail() const
{
    return solver_->inform();
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

const Eigen::VectorXd& LSSOLSolver::SI_result() const
{
    return solver_->result();
}

void LSSOLSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->resize(nrVar, nrEq + nrInEq, Eigen::lssol::QP2);
    Q_.resize(nrVar, nrVar);
    A_.resize(nrEq + nrInEq, nrVar);
    bl_.resize(nrEq + nrInEq);
    bu_.resize(nrEq + nrInEq);
}

bool LSSOLSolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
    const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
    const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    auto nrVar = Aeq.cols();
    auto nrEq = Aeq.rows();
    auto nrInEq = Aineq.rows();

    A_.block(0, 0, nrEq, nrVar) = Aeq;
    A_.block(nrEq, 0, nrInEq, nrVar) = Aineq;

    bl_.head(nrEq) = beq;
    bl_.tail(nrInEq).fill(-std::numeric_limits<double>::infinity());

    bu_.head(nrEq) = beq;
    bu_.tail(nrInEq) = bineq;

    Q_ = Q;

    return solver_->solve(XL, XU, Q_, c, A_, bl_, bu_);
}

} // namespace pc
