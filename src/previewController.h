#pragma once

#include "solverUtils.h"
#include <Eigen/Core>
#include <boost/timer/timer.hpp>
#include <memory>
#include <vector>

namespace pc
{

//Forward declaration
class Constrain;

/**
 * Structure representing all variables of a system for performing preview control.
 * Such system is defined as follow:
 * \f$X_{k+1} = Ax_{k} + Bu_{k} + d\f$.
 * After performing a recursion, this system can be represented in as:
 * \f$X_{k+1} = \Phi x_{0} + \Psi U + \Xi\f$
 */
struct PreviewSystemData
{
    PreviewSystemData(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                      const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                      int numberOfSteps);
    int nrStep; /**< The number of iteration to perform. */
    int xDim; /**< The dimension of the state vector */
    int uDim; /**< The dimension of the control vector */
    int fullXDim; /**< The full dimension of the state vector (xDim*nbStep) */
    int fullUDim; /**< The full dimension of the control vector (uDim*nbStep) */
    Eigen::VectorXd x0; /**< The initial state */
    Eigen::VectorXd xd; /**< The desired trajectory or desired final point */
    Eigen::MatrixXd A; /**< The state matrix */
    Eigen::MatrixXd B; /**< The control matrix */
    Eigen::VectorXd d; /**< The bias vector */
    Eigen::MatrixXd Phi; /**< The full (after nrStep) state matrix */
    Eigen::MatrixXd Psi; /**< The full (after nrStep) control matrix */
    Eigen::VectorXd xi; /**< The full (after nrStep) bias vector */
};

/**
 * The controller itself.
 * This class gives all the needed composants for performing a model preview control.
 * @warning This class waits for a discretized system ! Continuous systems are not yet implemented.
 */
class PreviewController
{
  public:
    /** 
     * Constructor of the class.
     * @param state The state matrix of the system.
     * @param control The control matrix of the system.
     * @param xInit The initial state.
     * @param xTraj The desired trajectory or final point.
     * @param numberOfSteps The number of step to perform.
     * @param pcFlag The type of model preview control @see pc::PCFlag.
     * @param pcFlag The solver to use @see pc::SolverFlag.
     */
    PreviewController(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                      const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                      int numberOfSteps, PCFlag pcFlag = PCFlag::Last, SolverFlag sFlag = SolverFlag::DEFAULT);
    /** 
     * Constructor of the class.
     * @param state The state matrix of the system.
     * @param control The control matrix of the system.
     * @param bias The bias matrix of the system.
     * @param xInit The initial state.
     * @param xTraj The desired trajectory or final point.
     * @param numberOfSteps The number of step to perform.
     * @param pcFlag The type of model preview control @see pc::PCFlag.
     * @param pcFlag The solver to use @see pc::SolverFlag.
     */
    PreviewController(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                      const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                      int numberOfSteps, PCFlag pcFlag = PCFlag::Last, SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
     * Select a solver. This function can be called at any time.
     * @param flag The solver to use @see pc::SolverFlag.
     */
    void selectQPSolver(SolverFlag flag);
    /**
     * Solve the system. 
     * @return True if a solution has been found.
     */
    bool solve();

    /**
     * Get the solver result.
     * @return The control vector \f$U\f$.
     */
    const Eigen::VectorXd &control() const noexcept;
    /**
     * Get the preview trajectory.
     * @return The trajectory vector \f$X\f$.
     */
    Eigen::VectorXd trajectory() const noexcept;
    /**
     * The time needed to solve the qp problem.
     * @return The elapsed time for solving.
     */
    boost::timer::cpu_times solveTime() const noexcept;
    /**
     * The time needed to build and solve the qp problem.
     * @return The elapsed time for build and solving.
     */
    boost::timer::cpu_times solveAndBuildTime() const noexcept;

    /**
     * Set the weights of the system.
     * All the weights will be the same.
     * @param wx Weight of the state.
     * @param wu Weight of the control.
     */
    void weights(double wx, double wu);
    /**
     * Set the weights of the system.
     * @param wx Weight of the state.
     * @param wu Weight of the control.
     */
    void weights(const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu);

    /**
     * Add a constrain to the system.
     * @param constr A constrain type @see TrajectoryConstrain. @see ControlConstrain.
     */
    void addConstrain(Constrain &constr);
    /**
     * Clear the constrains
     */
    void resetConstrains() noexcept;

  private:
    /**
     * Compute the preview system.
     */
    void previewSystem();
    /**
     * QP-like format.
     */
    void makeQPForm();

  private:
    PCFlag pcFlag_;
    int nrConstr_;

    std::unique_ptr<PreviewSystemData> psd_;
    std::vector<Constrain *> constr_;
    std::unique_ptr<SolverInterface> sol_;

    Eigen::MatrixXd Q_, AInEq_;
    Eigen::VectorXd c_, bInEq_, Wx_, Wu_;

    boost::timer::cpu_timer solveTime_, solveAndBuildTime_;
};

/**
 * Constrain to add to the system.
 * This is a pure virtual class
 */ 
class Constrain
{
  public:
    /**
     * Constructor of a constrain.
     * @param E The inequality constrain matrix.
     * @param f The inequality constrain vector.
     */
    Constrain(const Eigen::MatrixXd &E, const Eigen::VectorXd &f);

    /**
     * Initialization of the constrain.
     * @param psd An unique pointer to a PreviewSystemData.
     */
    virtual void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) = 0;
    /**
     * Update the constrain.
     * @param psd An unique pointer to a PreviewSystemData.
     */
    virtual void update(const std::unique_ptr<PreviewSystemData> &ps) = 0;

    /**
     * Function that return the number of constrains.
     * @return The number of constrains.
     */
    int nrConstr() noexcept
    {
        return nrConstr_;
    }
    /**
     * Function that return the QP form-like inequality matrix.
     * @return The inequality matrix.
     */
    const Eigen::MatrixXd &A() noexcept
    {
        return A_;
    }
    /**
     * Function that return the QP form-like inequality vector.
     * @return The inequality vector.
     */
    const Eigen::VectorXd &b() noexcept
    {
        return b_;
    }

  protected:
    int nrConstr_;
    Eigen::MatrixXd E_, A_;
    Eigen::VectorXd f_, b_;
};

class TrajectoryConstrain final : public Constrain
{
  public:
    /**
     * Constructor of a constrain.
     * Constrain of type \f$EX <= f\f$ 
     * @param E The inequality constrain matrix.
     * @param f The inequality constrain vector.
     */
    using Constrain::Constrain; //Inherits base class constructor
    void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) override;
    void update(const std::unique_ptr<PreviewSystemData> &psd) override;
};

class ControlConstrain final : public Constrain
{
  public:
    /**
     * Constructor of a constrain.
     * Constrain of type \f$EU <= f\f$ 
     * @param E The inequality constrain matrix.
     * @param f The inequality constrain vector.
     */
    using Constrain::Constrain; //Inherits base class constructor
    void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) override;
    void update(const std::unique_ptr<PreviewSystemData> &ps) override;
};

} // namespace pc