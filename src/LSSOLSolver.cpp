/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "LSSOLSolver.h"
#include <iostream>

namespace copra {

/*
 * LSSOL
 */

int LSSOLSolver::SI_fail() const
{
    return solver_.inform();
}

int LSSOLSolver::SI_iter() const
{
    return solver_.iter();
}

int LSSOLSolver::SI_maxIter() const
{
    return solver_.optimalityMaxIter();
}

void LSSOLSolver::SI_maxIter(int maxIter)
{
    solver_.optimalityMaxIter(maxIter);
    solver_.feasibilityMaxIter(maxIter);
}

void LSSOLSolver::SI_inform() const
{
    solver_.inform(std::cout);
}

void LSSOLSolver::SI_printLevel(int pl)
{
    solver_.printLevel(pl);
}

double LSSOLSolver::SI_feasibilityTolerance() const
{
    return solver_.feasibilityTol();
}

void LSSOLSolver::SI_feasibilityTolerance(double tol)
{
    solver_.feasibilityTol(tol);
}

bool LSSOLSolver::SI_warmStart() const
{
    return solver_.warm();
}

void LSSOLSolver::SI_warmStart(bool w)
{
    solver_.warm(w);
}

const Eigen::VectorXd& LSSOLSolver::SI_result() const
{
    return solver_.result();
}

void LSSOLSolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_.resize(nrVar, nrEq + nrInEq, Eigen::lssol::QP2);
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

    return solver_.solve(XL, XU, Q_, c, A_, bl_, bu_);
}

} // namespace copra
