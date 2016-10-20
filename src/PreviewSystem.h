#pragma once

#include <Eigen/Core>

namespace mpc
{

/**
 * Structure representing all variables of a system for performing preview control.
 * Such system is defined as follow:
 * \f$X_{k+1} = Ax_{k} + Bu_{k} + d\f$.
 * After performing a recursion, this system can be represented in as:
 * \f$X_{k+1} = \Phi x_{0} + \Psi U + \Xi\f$
 */
struct PreviewSystem
{
    /** 
     * Constructor of the class.
     * @param state The state matrix of the system.
     * @param control The control matrix of the system.
     * @param xInit The initial state.
     * @param xTraj The desired trajectory or final point.
     * @param numberOfSteps The number of step to perform.
     * @param sFlag The solver to use.
     */
    void system(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                int numberOfSteps);

    /** 
     * Constructor of the class.
     * @param state The state matrix of the system.
     * @param control The control matrix of the system.
     * @param xInit The initial state.
     * @param xTraj The desired trajectory or final point.
     * @param numberOfSteps The number of step to perform.
     * @param sFlag The solver to use.
     */
    void system(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                int numberOfSteps);

    int nrStep;          /**< The number of iteration to perform. */
    int xDim;            /**< The dimension of the state vector */
    int uDim;            /**< The dimension of the control vector */
    int fullXDim;        /**< The full dimension of the state vector (xDim*nbStep) */
    int fullUDim;        /**< The full dimension of the control vector (uDim*nbStep) */
    Eigen::VectorXd x0;  /**< The initial state */
    Eigen::VectorXd xd;  /**< The desired trajectory or desired final point */
    Eigen::MatrixXd A;   /**< The state matrix */
    Eigen::MatrixXd B;   /**< The control matrix */
    Eigen::VectorXd d;   /**< The bias vector */
    Eigen::MatrixXd Phi; /**< The full (after nrStep) state matrix */
    Eigen::MatrixXd Psi; /**< The full (after nrStep) control matrix */
    Eigen::VectorXd xi;  /**< The full (after nrStep) bias vector */
};

} // namespace mpc