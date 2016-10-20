#include "QLDSolver.h"

namespace mpc
{

/*
 * QLD
 */

QLDSolver::QLDSolver()
    : solver_(std::make_unique<Eigen::QLD>())
{
}

int QLDSolver::SI_fail() const
{
    return solver_->fail();
}

void QLDSolver::SI_inform() const
{
    solver_->inform();
}

void QLDSolver::SI_printLevel(int pl) const
{
    solver_->verbose(static_cast<bool>(pl));
}

void QLDSolver::SI_tol(double tol) const
{
    solver_->accuracyTol(tol);
}

const Eigen::VectorXd &QLDSolver::SI_result() const
{
    return solver_->result();
}

void QLDSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool QLDSolver::SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
                         const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
                         const Eigen::MatrixXd &Aineq,
                         const Eigen::VectorXd &Bineq,
                         const Eigen::VectorXd &XL, const Eigen::VectorXd &XU)
{
    return solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq, XL, XU);
}

} // namespace pc