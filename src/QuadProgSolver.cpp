#include "QuadProgSolver.h"

namespace pc
{

/*
 * QuadProg dense
 */

QuadProgDenseSolver::QuadProgDenseSolver()
    : solver_(std::make_unique<Eigen::QuadProgDense>())
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
    //QuadProg does not have bound constrains.
    //They will be transformed into inequality constrains.
    //The bound limits constrain size is 2 times (upper and lower bound) the number of variables
    solver_->problem(nrVar, nrEq, nrInEq + 2 * nrVar);
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

} // namespace pc