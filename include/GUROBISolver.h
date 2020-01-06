/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "SolverInterface.h"
#include <Eigen/Core>
#include <eigen-gurobi/Gurobi.h>

namespace copra {

/**
 * GUROBISolver solver for both dense matrix.
 */
class COPRA_DLLAPI GUROBISolver : public SolverInterface // TODO: Enable sparse matrix
{
public:
    /**
     * GUROBISolver default constructor
     */
    GUROBISolver();

    int SI_fail() const override;
    void SI_inform() const override;
    int SI_iter() const override;
    void SI_printLevel(int pl) override;
    void SI_feasibilityTolerance(double tol) override;
    int SI_maxIter() const override;
    void SI_maxIter(int maxIter) override;

    bool SI_warmStart() const override;
    void SI_warmStart(bool w) override;

    /**
     * Get the solver's solution.
     * \return The qp solver result.
     */
    const Eigen::VectorXd& SI_result() const override;

    /**
     * Initialize the variables of the problem to solve.
     * \see SolverInterface::SI_problem()
     * \return The qp solver result.
     */
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;

    /**
     * Solve the problem.
     * \see SolverInterface::SI_solve()
     * \return The qp solver result.
     */
    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;

    Eigen::GurobiDense& baseSolver() noexcept { return solver_; }

private:
    Eigen::GurobiDense solver_;
};

} // namespace pc
