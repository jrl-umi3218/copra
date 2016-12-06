// This file is part of ModelPreviewController.

// ModelPreviewController is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ModelPreviewController is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ModelPreviewController.  If not, see
// <http://www.gnu.org/licenses/>.

#pragma once

#include "solverUtils.h"
#include <Eigen/Core>
#include <boost/timer/timer.hpp>
#include <memory>
#include <vector>
#include <string>

namespace mpc
{

// Forward declaration
class Constraint;
class EqIneqConstraint;
class ControlBoundConstraint;
struct PreviewSystem;
enum class ConstraintFlag;

/**
 * The controller itself.
 * This class gives all the needed composants for performing a model preview
 * control.
 * @warning This class waits for a discretized system ! Continuous systems are
 * not yet implemented.
 */
class MPCTypeFull
{
  public:
    /**
   * Initialize problem variables to default and get the desired solver
   * You need to call initializeController before using the MPCTypeFull
   * @param ps A preview system to amke a copy from.
   * @param sFlag The flag corresponding to the desired solver.
   */
    MPCTypeFull(SolverFlag sFlag = SolverFlag::DEFAULT);
    /**
   * Initialize problem variables w.r.t. the PreviewSystem and get the desired
   * solver
   * @param ps A preview system to amke a copy from.
   * @param sFlag The flag corresponding to the desired solver.
   */
    MPCTypeFull(const std::shared_ptr<PreviewSystem> &ps,
                SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
   * Select a solver. This function can be called at any time.
   * @param flag The solver to use @see pc::SolverFlag.
   */
    void selectQPSolver(SolverFlag flag);

    /**
   * Initialize the controller w.r.t. the preview system.
   * This function needs to be called each time the system dimension changes.
   * @param ps The preview system
   */
    virtual void initializeController(const std::shared_ptr<PreviewSystem> &ps);

    /**
   * Solve the system.
   * @return True if a solution has been found.
   * Fill Phi, Psi, xi in PreviewSystem
   * Fill A, b in Constrains
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
   * @param wx Weight of the state.
   * @param wu Weight of the control.
   */
    virtual void weights(const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu);

    /**
   * Add a constraint to the system. The shared_ptr if not copied !
   * So, if it deleted before solving, the MPC will not use the constraint.
   * In case of a unwilling deletion, a warning is displayed in Debug mode.
   * @param constr A constraint type @see TrajectoryConstrain. @see ControlConstrain.
   * @see TrajectoryBoundConstraint @see ControlBoundConstraint
   */
    void addConstraint(const std::shared_ptr<Constraint>& constr);

    /**
     * Clear the constrains
     */
    void resetConstraints() noexcept;

  protected:
    /**
     */
    void addConstraintByType(const std::shared_ptr<Constraint>& constr);

    /**
     * Update the system and its constrains.
     * @param ps The preview system
     * Fill Phi, Psi, xi in PreviewSystem
     * Fill A, b in Constrains
     */
    void updateSystem();

    /**
     * QP-like format.
     */
    virtual void makeQPForm();

    /**
     * Check if the constraints still exist.
     * Output into std::cerr if something wrong happened.
     */
    void checkConstraints();

  protected:
    struct Constraints
    {
        Constraints();
        void clear();
        void updateNr();

        int nrConstr;
        int nrEqConstr;
        int nrIneqConstr;
        int nrBoundConstr;
        std::vector<std::pair<std::weak_ptr<Constraint>, std::string>> wpConstr;
        std::vector<std::pair<std::weak_ptr<EqIneqConstraint>, std::string>> wpEqConstr;
        std::vector<std::pair<std::weak_ptr<EqIneqConstraint>, std::string>> wpIneqConstr;
        std::vector<std::pair<std::weak_ptr<ControlBoundConstraint>, std::string>> wpBoundConstr;
    };

  protected:
    std::shared_ptr<PreviewSystem> ps_;
    std::unique_ptr<SolverInterface> sol_;
    Constraints constraints_;

    Eigen::MatrixXd Q_, Aineq_, Aeq_;
    Eigen::VectorXd c_, bineq_, beq_, lb_, ub_, Wx_, Wu_;

    boost::timer::cpu_timer solveTime_, solveAndBuildTime_;
};

/**
 * This is a variant of MPCTypeFull controller
 * In spite of using the full generated preview matrices, it uses only the last
 * part (the most in the future part).
 * Thus, this class has way faster results with the disadvantages of not
 * considering a minimization along all the trajectory.
 * @warning This class waits for a discretized system ! Continuous systems are
 * not yet implemented.
 */
class MPCTypeLast final : public MPCTypeFull
{
  public:
    /**
   * See @see MPCTypeFull::MPCTypeFull
   */
    MPCTypeLast(SolverFlag sFlag = SolverFlag::DEFAULT);
    /**
   * See @see MPCTypeFull::MPCTypeFull
   */
    MPCTypeLast(const std::shared_ptr<PreviewSystem> &ps, SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
   * See @see  MPCTypeFull::initializeController
   */
    void initializeController(const std::shared_ptr<PreviewSystem> &ps) override;

    /**
   * See @see MPCTypeFull::weights
   */
    void weights(const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu) override;

  protected:
    void makeQPForm() override;
};

} // namespace pc