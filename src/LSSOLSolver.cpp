#include "LSSOLSolver.h"

namespace pc {

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

void LSSOLSolver::SI_printLevel(int pl) const
{
    solver_->printLevel(pl);
}

void LSSOLSolver::SI_tol(double tol) const
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
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool LSSOLSolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& C,
    const Eigen::MatrixXd& Aeq,
    const Eigen::VectorXd& Beq,
    const Eigen::MatrixXd& Aineq,
    const Eigen::VectorXd& Bineq,
    const Eigen::VectorXd& XL,
    const Eigen::VectorXd& XU)
{
    return solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq, XL, XU);
}

} // namespace pc