#pragma once

#include "solverUtils.h"
#include <Eigen/Core>
#include <boost/timer/timer.hpp>
#include <memory>
#include <vector>

namespace mpc
{

//Forward declaration
class Constrain;
class PreviewSystem;

/**
 * The controller itself.
 * This class gives all the needed composants for performing a model preview control.
 * @warning This class waits for a discretized system ! Continuous systems are not yet implemented.
 */
class MPCTypeFull
{
  public:
    /**
     * Make a copy of the PreviewSystem and get the desired solver
     * @param ps A preview system to amke a copy from.
     * @param sFlag The flag corresponding to the desired solver.
     */
    MPCTypeFull(const PreviewSystem &ps, SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
     * Select a solver. This function can be called at any time.
     * @param flag The solver to use @see pc::SolverFlag.
     */
    void selectQPSolver(SolverFlag flag);

    /**
	 * Update the system and its constrains.
	 * @param ps The preview system
	 * Fill Phi, Psi, xi in PreviewSystem
	 * Fill A, b in Constrains
	 */
    void updateSystem(PreviewSystem &ps);

    /**
     * Solve the system. 
     * @return True if a solution has been found.
     */
    bool solve(const PreviewSystem &ps);

    /**
     * Get the solver result.
     * @return The control vector \f$U\f$.
     */
    const Eigen::VectorXd &control() const noexcept;
    /**
     * Get the preview trajectory.
     * @return The trajectory vector \f$X\f$.
     */
    Eigen::VectorXd trajectory(const PreviewSystem &ps) const noexcept;
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
     * @param wx Weight of the state.
     * @param wu Weight of the control.
     */
    virtual void weights(const PreviewSystem &ps, const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu);

    /**
     * Add a constrain to the system.
     * @param constr A constrain type @see TrajectoryConstrain. @see ControlConstrain.
     */
    void addConstrain(const PreviewSystem &ps, Constrain &constr);
    /**
     * Clear the constrains
     */
    void resetConstrains() noexcept;

  protected:
    /**
     * Compute the preview system.
     */
    void previewSystem(const PreviewSystem &ps);
    /**
     * QP-like format.
     */
    virtual void makeQPForm(const PreviewSystem &ps);

  protected:
    int nrConstr_;

    std::vector<Constrain *> constr_;
    std::unique_ptr<SolverInterface> sol_;

    Eigen::MatrixXd Q_, AInEq_;
    Eigen::VectorXd c_, bInEq_, Wx_, Wu_;

    boost::timer::cpu_timer solveTime_, solveAndBuildTime_;
};

/**
 * This is a variant of MPCTypeFull controller
 * In spite of using the full generated preview matrices, it uses only the last part (the most in the future part).
 * Thus, this class has way faster results with the disadvantages of not considering a minimization along all the trajectory.
 * @warning This class waits for a discretized system ! Continuous systems are not yet implemented.
 */
class MPCTypeLast final : public MPCTypeFull
{
  public:
    /**
     * See @see MPCTypeFull::MPCTypeFull
     */
    MPCTypeLast(const PreviewSystem &ps, SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
     * See @see MPCTypeFull::weights
     */
    void weights(const PreviewSystem &ps, const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu) override;

  protected:
    void makeQPForm(const PreviewSystem &ps) override;
};

} // namespace pc