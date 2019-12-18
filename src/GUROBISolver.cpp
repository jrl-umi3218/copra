/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "GUROBISolver.h"
#include <iostream>

namespace copra {

/*
 * GUROBI
 */

GUROBISolver::GUROBISolver()
{
    solver_.displayOutput(false);
    solver_.feasibilityTolerance(1e-8);
}

int GUROBISolver::SI_fail() const
{
    return solver_.fail();
}

int GUROBISolver::SI_iter() const
{
    return solver_.iter();
}

void GUROBISolver::SI_printLevel(int pl)
{
    solver_.displayOutput(pl != 0);
}

double GUROBISolver::SI_feasibilityTolerance() const
{
    solver_.feasibilityTolerance();
}

void GUROBISolver::SI_feasibilityTolerance(double tol)
{
    solver_.feasibilityTolerance(tol); //primal feasible tolerance
}

int GUROBISolver::SI_maxIter() const
{
    solver_.optimalityTolerance();
}

void GUROBISolver::SI_maxIter(int maxIter)
{
    solver_.optimalityTolerance(tol); //dual feasible tolerance
}

bool GUROBISolver::SI_warmStart() const
{
    // 2 is for non-warmstart solver
    return solver_.warmStart() != Eigen::GurobiCommon::WarmStatus::NONE;
}

void GUROBISolver::SI_warmStart(bool w)
{
    // -1 is is for default warmstart
    // 2 is for non-warmstart solver
    using WS = Eigen::GurobiCommon::WarmStatus;
    solver_.warmStart((w ? WS::DEFAULT : WS::NONE));
}

void GUROBISolver::SI_inform() const
{
    solver_.inform();
}

const Eigen::VectorXd& GUROBISolver::SI_result() const
{
    return solver_.result();
}

void GUROBISolver::SI_problem(int nrVar, int nrEq, int nrInEq)
{
    solver_.problem(nrVar, nrEq, nrInEq);
}

bool GUROBISolver::SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
    const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
    const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
    const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
{
    return solver_.solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU);
}

} // namespace pc
