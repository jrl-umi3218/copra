#include "solvers.h"

// std
#include <utility>

namespace pc
{

/*
 * SolverInterface
 */
int SolverInterface::SI_fail() const
{
    return 0;
}

const Eigen::VectorXd &SolverInterface::SI_result() const
{
    return std::move(Eigen::VectorXd());
}

void SolverInterface::SI_problem(int /* nrVar */, int /* nrEq */, int /* nrInEq */)
{
}

bool SolverInterface::SI_solve(const Eigen::MatrixXd & /* Q */, const Eigen::VectorXd & /* C */,
                               const Eigen::MatrixXd & /* Aeq */, const Eigen::VectorXd & /* Beq */,
                               const Eigen::MatrixXd & /* Aineq */, const Eigen::VectorXd & /* Bineq */,
                               const Eigen::VectorXd & /* XL */, const Eigen::VectorXd & /* XU */)
{
    return false;
}

/*
 * QuadProg dense
 */

QuadProgDenseSolver::QuadProgDenseSolver()
    : solver_(std::make_unique<Eigen::QuadProgDense>())
{
}

QuadProgDenseSolver::QuadProgDenseSolver(int nrVar, int nrEq, int nrInEq)
    : solver_(std::make_unique<Eigen::QuadProgDense>(nrVar, nrEq, nrInEq))
{
}

int QuadProgDenseSolver::SI_fail() const
{
    return solver_->fail();
}

const Eigen::VectorXd &QuadProgDenseSolver::SI_result() const
{
    return solver_->result();
}

void QuadProgDenseSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_->problem(nrVar, nrEq, nrInEq);
}

bool QuadProgDenseSolver::SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
                                   const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
                                   const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq,
                                   const Eigen::VectorXd &XL, const Eigen::VectorXd &XU)
{
    auto nrLines = XL.rows();
    auto ALines = Aineq.rows();
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nrLines, nrLines);
    Eigen::MatrixXd ineqMat(ALines + 2 * nrLines, Aineq.cols());
    Eigen::VectorXd ineqVec(ALines + 2 * nrLines);
    ineqMat.topRows(ALines) = Aineq;
    ineqMat.block(ALines, 0, nrLines, nrLines) = I;
    ineqMat.bottomRows(nrLines) = -I;
    ineqVec.head(ALines) = Bineq;
    ineqVec.segment(ALines, nrLines) = XU;
    ineqVec.tail(nrLines) = -XL;

    return solver_->solve(Q, C, Aeq, Beq, ineqMat, ineqVec);
}

/*
 * QLD
 */

QLDSolver::QLDSolver()
    : solver_(std::make_unique<Eigen::QLD>())
{
}

QLDSolver::QLDSolver(int nrVar, int nrEq, int nrInEq)
    : solver_(std::make_unique<Eigen::QLD>(nrVar, nrEq, nrInEq))
{
}

int QLDSolver::SI_fail() const
{
    return solver_->fail();
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
                         const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq,
                         const Eigen::VectorXd &XL, const Eigen::VectorXd &XU)
{
    return solver_->solve(Q, C, Aeq, Beq, Aineq, Bineq, XL, XU, 1e-3); //TODO: Change harcoded value
}

} // namespace pc