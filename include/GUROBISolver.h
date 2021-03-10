/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "SolverInterface.h"
#include <Eigen/Core>
#include <eigen-gurobi/Gurobi.h>

namespace copra {

/*! \brief GUROBISolver solver for both dense matrix.*/
class COPRA_DLLAPI GUROBISolver : public SolverInterface { // TODO: Enable sparse matrix
public:
    /*! \brief Default constructor. */
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
    const Eigen::VectorXd& SI_result() const override;
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;
    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;
    /*! \brief Return underlying solver. */
    Eigen::GurobiDense& baseSolver() noexcept { return solver_; }

private:
    Eigen::GurobiDense solver_; /*!< QP solver */
};

} // namespace copra
