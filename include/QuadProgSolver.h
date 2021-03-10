/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "SolverInterface.h"
#include <Eigen/Core>
#include <eigen-quadprog/QuadProg.h>

namespace copra {

/*! \brief QuadProg solver for dense matrix. */
class COPRA_DLLAPI QuadProgDenseSolver : public SolverInterface {
public:
    /*! \brief Default constructor */
    QuadProgDenseSolver() = default;

    /*! \brief Get information of eventual fail's solver output as define by the
     * solver documentation.
     * \return 0 No problems
     * \return 1 The minimization problem has no solution
     * \return 2 Problems with the decomposition of Q (Is it symmetric positive
     * definite matrix?)
     */
    int SI_fail() const override;
    int SI_iter() const override;
    void SI_inform() const override;
    const Eigen::VectorXd& SI_result() const override;
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;
    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;
    /*! \brief Return underlying solver. */
    Eigen::QuadProgDense& baseSolver() noexcept { return solver_; }

private:
    Eigen::QuadProgDense solver_; /*!< QP solver */
};

} // namespace copra
