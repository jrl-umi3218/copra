// This file is part of mpc.

// mpc is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with mpc.  If not, see
// <http://www.gnu.org/licenses/>.

#pragma once

// stl
#include <memory>
#include <string>
#include <vector>

// boost
#include <boost/timer/timer.hpp>

// eigen
#include <Eigen/Core>

// mpc
#include "config.hh"
#include "debugUtils.h"
#include "solverUtils.h"
#include "typedefs.h"

namespace mpc {

// Forward declaration
class Constraint;
class EqIneqConstraint;
class ControlBoundConstraint;
struct PreviewSystem;
enum class ConstraintFlag;

/**
 * The controller itself.
 * This class gives all the needed composants for performing a model preview control.
 * It solves:\n
 * \f$X = \Phi x_{0} + \Psi U + \Xi\f$, where \f$U\f$ is the optimization vector.
 * \note \f$X = [x_0^T x_1^T ... x_N^T]^T\f$ and \f$U = [u_0^T u_1^T ... u_{N-1}^T]^T\f$
 * where \f$N\f$ is the dimension of the system (the number of steps).
 * \warning This class waits for a discretized system ! Continuous systems are not implemented.
 */
class MPC_DLLAPI MPCTypeFull {
public:
    /**
     * Initialize problem variables to default and get the desired solver
     * You need to call initializeController before using the MPCTypeFull
     * \param ps A preview system to amke a copy from.
     * \param sFlag The flag corresponding to the desired solver.
     */
    MPCTypeFull(SolverFlag sFlag = SolverFlag::DEFAULT);
    /**
     * Initialize problem variables w.r.t. the PreviewSystem and get the desired
     * solver
     * \param ps A preview system to amke a copy from.
     * \param sFlag The flag corresponding to the desired solver.
     */
    MPCTypeFull(const std::shared_ptr<PreviewSystem>& ps,
        SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
     * Select a solver. This function can be called at any time.
     * \param flag The solver to use \see pc::SolverFlag.
     */
    void selectQPSolver(SolverFlag flag);

    /**
     * Initialize the controller with regard to the preview system.
     * This function needs to be called each time the system dimension changes.
     * \param ps The preview system
     */
    virtual void initializeController(const std::shared_ptr<PreviewSystem>& ps);

    /**
     * Solve the system.
     * \return True if a solution has been found.
     * Fill Phi, Psi, xi in PreviewSystem
     * Fill A, b in Constraints
     */
    bool solve();

    /**
     * Get the solver result.
     * \return The control vector \f$U\f$.
     */
    const Eigen::VectorXd& control() const noexcept;
    /**
     * Get the preview trajectory.
     * \return The trajectory vector \f$X\f$.
     */
    Eigen::VectorXd trajectory() const noexcept;
    /**
     * The time needed to solve the qp problem.
     * \return The elapsed time for solving.
     */
    boost::timer::cpu_times solveTime() const noexcept;
    /**
     * The time needed to build and solve the qp problem.
     * \return The elapsed time for build and solving.
     */
    boost::timer::cpu_times solveAndBuildTime() const noexcept;

    /**
     * Set the weights of the system.
     * Perform a move semantic if a fullsize vector is given as rvalue (this is faster).
     * \param wx Weight of the state.
     * \param wu Weight of the control.
     * \throw Throw a std::domain_error is Wx or Wu are badly dimension.
     */
    template <typename TVec1, typename TVec2,
        typename = std::enable_if_t<!is_all_arithmetic<TVec1, TVec2>::value> >
    void weights(TVec1&& Wx, TVec2&& Wu)
    {
        if (Wx.rows() == Wx_.rows())
            Wx_ = std::forward<TVec1>(Wx);
        else if (Wx.rows() == ps_->xDim)
            for (auto i = 0; i < ps_->nrXStep; ++i)
                Wx_.segment(i * ps_->xDim, ps_->xDim) = Wx;
        else
            DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSXDim("Wx", Wx, ps_.get()));

        if (Wu.rows() == ps_->fullUDim)
            Wu_ = std::forward<TVec2>(Wu);
        else if (Wu.rows() == ps_->uDim)
            for (auto i = 0; i < ps_->nrUStep; ++i)
                Wu_.segment(i * ps_->uDim, ps_->uDim) = Wu;
        else
            DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSUDim("Wu", Wu, ps_.get()));
    }

    /**
     * Set the weights of the system. All variables are set to the same weight.
     * \param wx Weight of the state.
     * \param wu Weight of the control.
     */
    template <typename T1, typename T2,
        typename = std::enable_if_t<is_all_arithmetic<T1, T2>::value> >
    void weights(T1 Wx, T2 Wu)
    {
        Wx_.setConstant(Wx);
        Wu_.setConstant(Wu);
    }

    /**
     * Add a constraint to the system. The shared_ptr if not copied !
     * So, if it deleted before solving, the MPC will not use the constraint.
     * In case of a unwilling deletion, a warning is displayed in Debug mode.
     * \param constr A constraint type \see TrajectoryConstrain. \see ControlConstrain.
     * \see TrajectoryBoundConstraint \see ControlBoundConstraint
     */
    void addConstraint(const std::shared_ptr<Constraint>& constr);

    /**
     * Clear the constraints
     */
    void resetConstraints() noexcept;

protected:
    /**
     * Add constraints into constraints_ \see Constraints
     * \param constr The constraint to add
     */
    void addConstraintByType(const std::shared_ptr<Constraint>& constr);

    /**
     * Resize Aeq, beq, Aineq, bineq, ub, lb to default.
     */
    void clearConstraintMatrices();

    /**
     * Update the system and its constraints.
     * \param ps The preview system
     * Fill A, b in Constraints
     */
    void updateSystem();

    /**
     * QP-like format.
     */
    virtual void makeQPForm();

    /**
     * Check if the constraints still exist.
     * In Debug mode: Output into std::cerr if a constraint has been deleted.
     * In Release mode: No output.
     */
    void checkDeleteConstraints();

protected:
    struct Constraints {
        Constraints();
        void clear();
        void updateNr();

        int nrEqConstr;
        int nrIneqConstr;
        std::vector<std::shared_ptr<Constraint> > spConstr;
        std::vector<std::shared_ptr<EqIneqConstraint> > spEqConstr;
        std::vector<std::shared_ptr<EqIneqConstraint> > spIneqConstr;
        std::vector<std::shared_ptr<ControlBoundConstraint> > spBoundConstr;
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
 * \warning This class waits for a discretized system ! Continuous systems are
 * not implemented.
 */
class MPC_DLLAPI MPCTypeLast : public MPCTypeFull {
public:
    MPCTypeLast(SolverFlag sFlag = SolverFlag::DEFAULT);
    MPCTypeLast(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag = SolverFlag::DEFAULT);

    void initializeController(const std::shared_ptr<PreviewSystem>& ps) override;

protected:
    void makeQPForm() override;
};

} // namespace pc