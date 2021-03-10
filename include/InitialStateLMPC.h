/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "LMPC.h"
#include <Eigen/LU>

namespace copra {

/*! \brief LMPC with optimization of initial state.
 * This class gives all the needed composants for performing a model preview control including optimization of the initial state.
 * \warning This class waits for a discretized system ! Continuous systems are not implemented.
 */
class COPRA_DLLAPI InitialStateLMPC : public LMPC {
public:
    InitialStateLMPC(SolverFlag sFlag = SolverFlag::DEFAULT);
    InitialStateLMPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag = SolverFlag::DEFAULT);

    /*! \brief The initial-state vector \f$x_{0}\f$. */
    Eigen::VectorXd initialState() const noexcept;
    /*! \brief Add a cost related to the initial state. */
    void resetInitialStateCost(const Eigen::MatrixXd& R, const Eigen::VectorXd& r);
    /*! \brief Add lower and upper bounds for the initial state. */
    void resetInitialStateBounds(const Eigen::VectorXd& l, const Eigen::VectorXd& u);

private:
    void clearConstraintMatrices() override;
    void updateQPMatrixSize() override;
    void updateConstraintMatrixSize() override;
    void makeQPForm() override;
    void updateResults() override;

protected:
    Eigen::MatrixXd R_; /*!< Initial state cost matrix. Must be a positive finite matrix. Default to \f$I \cdot 1e^{\minus 6}\f$ */
    Eigen::VectorXd r_; /*!< Initial state cost vector */
    Eigen::VectorXd x0lb_; /*!< Initial state lower bound */
    Eigen::VectorXd x0ub_; /*!< Initial state upper bound */
};

} // namespace copra
