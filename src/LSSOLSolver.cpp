#include "LSSOLSolver.h"

namespace pc
{

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

const Eigen::VectorXd &LSSOLSolver::SI_result() const
{
    return solver_->result();
}

void LSSOLSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool LSSOLSolver::SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
                           const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
                           const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq,
                           const Eigen::VectorXd &XL, const Eigen::VectorXd &XU)
{
    return solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq, XL, XU);
}

} // namespace pc