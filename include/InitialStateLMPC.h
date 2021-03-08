/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "LMPC.h"
#include <Eigen/LU>

namespace copra {


/**
 * The controller itself.
 * This class gives all the needed composants for performing a model preview control including optimization of the initial state.
 * \warning This class waits for a discretized system ! Continuous systems are not implemented.
 */
class COPRA_DLLAPI InitialStateLMPC : public LMPC {
public:
    /**
     * Initialize problem variables to default and get the desired solver
     * You need to call initializeController before using the MPCTypeFull
     * \param sFlag The flag corresponding to the desired solver.
     */
    InitialStateLMPC(SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
     * Initialize problem variables w.r.t. the PreviewSystem and get the desired
     * solver
     * \param ps A preview system to make a copy from.
     * \param sFlag The flag corresponding to the desired solver.
     */
    InitialStateLMPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
     * Initialize the controller with regard to the preview system.
     * This function needs to be called each time the system dimension changes.
     * \param ps The preview system
     */
    void initializeController(const std::shared_ptr<PreviewSystem>& ps) override;

    /**
     * Solve the system.
     * \return True if a solution has been found.
     * Fill Phi, Psi, xi in PreviewSystem
     * Fill A, b in Constraints
     */
    bool solve() override;

    /**
     * Get the solver result.
     * \return The initial-state vector \f$x_{0}\f$.
     */
    Eigen::VectorXd initialState() const noexcept;
    
    /**
     * Get the solver result.
     * \return The control vector \f$U\f$.
     */
    const Eigen::VectorXd& control() const noexcept override;

    /**
     * Get the preview trajectory.
     * \return The trajectory vector \f$X\f$.
     */
    Eigen::VectorXd trajectory() const noexcept override;

    /**
     * Add a cost related to the initial state
     */
    void resetInitialStateCost(const Eigen::MatrixXd& R, const Eigen::VectorXd& r);

    /**
     * Add lower and upper bounds for the initial state 
     */
    void resetInitialStateBounds(const Eigen::VectorXd& l, const Eigen::VectorXd& u);

protected:

    /**
     * Resize internal matrices and vectors to default.
     */
    void clearConstraintMatrices() override;

    /**
     * Update the system and its constraints.
     */
    void updateSystem() override;

    /**
     * QP-like format.
     */
    void makeQPForm() override;

protected:
    Eigen::MatrixXd R_;
    Eigen::VectorXd r_;
    Eigen::MatrixXd E_;
    Eigen::VectorXd f_;
    Eigen::VectorXd l_; //initial state lower bound
    Eigen::VectorXd u_; //initial state upper bound

    Eigen::MatrixXd Yeq_;
    Eigen::VectorXd zeq_;
    Eigen::MatrixXd Yineq_;
    Eigen::VectorXd zineq_;

    Eigen::MatrixXd newQ_;
    Eigen::VectorXd newc_;
    Eigen::VectorXd newlb_;
    Eigen::VectorXd newub_;

    Eigen::MatrixXd newAeq_;
    Eigen::VectorXd newbeq_;
    Eigen::MatrixXd newAineq_;
    Eigen::VectorXd newbineq_;

    Eigen::VectorXd control_;
};

} // namespace pc
