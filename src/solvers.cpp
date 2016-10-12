#include "solvers.h"

// std
#include <utility>

namespace pc
{

/*
 * QuadProg dense
 */

QuadProgDenseSolver::QuadProgDenseSolver()
    : solver_(std::make_unique<QuadProgDense>())
{
}

QuadProgDenseSolver::QuadProgDenseSolver(int nrVar, int nrEq, int nrInEq)
    : solver_(std::make_unique<QuadProgDense>(nrVar, nrEq, nrInEq))
{
}

const Eigen::VectorXi &QuadProgDenseSolver::SI_iter() const
{
    return solver_->iter();
}

int QuadProgDenseSolver::SI_fail() const
{
    return solver_->fail();
}

const Eigen::VectorXd &QuadProgDenseSolver::SI_result() const
{
    return seolver_->result();
}

void QuadProgDenseSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool QuadProgDenseSolver::SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
                                   const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
                                   const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq)
{
    solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq);
}

/*
 * QuadProg sparse
 */

QuadProgSparseSolver::QuadProgSparseSolver()
    : solver_(std::make_unique<QuadProgSparse>())
{
}

QuadProgSparseSolver::QuadProgSparseSolver(int nrVar, int nrEq, int nrInEq)
    : solver_(std::make_unique<QuadProgSparse>(nrVar, nrEq, nrInEq))
{
}

const Eigen::VectorXi &QuadProgSparseSolver::SI_iter() const
{
    return solver_->iter();
}

int QuadProgSparseSolver::SI_fail() const
{
    return solver_->fail();
}

const Eigen::VectorXd &QuadProgSparseSolver::SI_result() const
{
    return seolver_->result();
}

void QuadProgSparseSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool QuadProgSparseSolver::SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
                                    const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
                                    const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq)
{
    solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq);
}

/*
 * QLD
 */

QLDSolver::QLDSolver()
    : solver_(std::make_unique<QuadProgSparse>())
{
}

QLDSolver::QLDSolver(int nrVar, int nrEq, int nrInEq)
    : solver_(std::make_unique<QuadProgSparse>(nrVar, nrEq, nrInEq))
{
}

const Eigen::VectorXi &QLDSolver::SI_iter() const
{
    return solver_->iter();
}

int QLDSolver::SI_fail() const
{
    return solver_->fail();
}

const Eigen::VectorXd &QLDSolver::SI_result() const
{
    return seolver_->result();
}

void QLDSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool QLDSolver::SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
                         const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
                         const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq)
{
    solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq);
}

} // namespace pc