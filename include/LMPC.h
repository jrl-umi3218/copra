/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "debugUtils.h"
#include "solverUtils.h"
#include "typedefs.h"
#include <Eigen/Core>
#include <chrono>
#include <memory>
#include <string>
#include <vector>

namespace copra {

// Forward declaration
enum class ConstraintFlag;
struct PreviewSystem;
class Constraint;
class ControlBoundConstraint;
class EqIneqConstraint;
class CostFunction;

/**
 * The controller itself.
 * This class gives all the needed composants for performing a model preview control.
 * It solves:\n
 * \f$X = \Phi x_{0} + \Psi U + \Xi\f$, where \f$U\f$ is the optimization vector.
 * \note \f$X = [x_0^T x_1^T ... x_N^T]^T\f$ and \f$U = [u_0^T u_1^T ... u_{N-1}^T]^T\f$
 * where \f$N\f$ is the dimension of the system (the number of steps).
 * \warning This class waits for a discretized system ! Continuous systems are not implemented.
 */
class COPRA_DLLAPI LMPC {
public:
    /**
     * Initialize problem variables to default and get the desired solver
     * You need to call initializeController before using the MPCTypeFull
     * \param sFlag The flag corresponding to the desired solver.
     */
    LMPC(SolverFlag sFlag = SolverFlag::DEFAULT);
    /**
     * Initialize problem variables w.r.t. the PreviewSystem and get the desired
     * solver
     * \param ps A preview system to amke a copy from.
     * \param sFlag The flag corresponding to the desired solver.
     */
    LMPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag = SolverFlag::DEFAULT);

    /**
     * Select a solver. It will load a solver with default values.
     * \see useSolver if you want to give to the MPC a solver with specific parameters
     * \param flag The solver to use \see pc::SolverFlag.
     */
    void selectQPSolver(SolverFlag flag);

    /**
     * Make the MPC uses a user-defined qp solver.
     * \param solver The user-defined solver. It must inheritate from the SolverInterface.
     */
    void useSolver(std::unique_ptr<SolverInterface>&& solver);

    /**
     * Initialize the controller with regard to the preview system.
     * This function needs to be called each time the system dimension changes.
     * \param ps The preview system
     */
    void initializeController(const std::shared_ptr<PreviewSystem>& ps);

    /**
     * Solve the system.
     * \return True if a solution has been found.
     * Fill Phi, Psi, xi in PreviewSystem
     * Fill A, b in Constraints
     */
    bool solve();

    /**
     * Print information on the QP solver status.
     */
    void inform() const noexcept;

    /**
     * Get a reference to the solver if any.
     * \return QP solver currently used.
     */
    inline SolverInterface& solver() { return *sol_; }
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
     * \return The elapsed time (in s) for solving.
     */
    double solveTime() const noexcept;
    /**
     * The time needed to build and solve the qp problem.
     * \return The elapsed time (in s) for build and solving.
     */
    double solveAndBuildTime() const noexcept;

    /**
     * Add a cost function to the system.
     * \param costFun A cost type \see TrajectoryCost \see TargetCost
     * \see ControlCost \see MixedTrajectoryCost \see MixedTargetCost
     */
    void addCost(const std::shared_ptr<CostFunction>& costFun);

    /**
     * Add a constraint to the system.
     * \param constr A constraint type \see TrajectoryConstrain \see ControlConstraint
     * \see TrajectoryBoundConstraint \see ControlBoundConstraint
     */
    void addConstraint(const std::shared_ptr<Constraint>& constr);

    /**
     * Clear the constraints
     */
    void resetConstraints() noexcept;

    /**
     * Remove cost
     */
    void removeCost(const std::shared_ptr<CostFunction>& costFun);

    /**
     * Remove cost
     */
    void removeConstraint(const std::shared_ptr<Constraint>& constr);

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
     * Fill A, b in Constraints
     */
    void updateSystem();

    /**
     * QP-like format.
     */
    void makeQPForm();

    /**
     * Check if a cost or a constraint still exist.
     * In Debug mode: Output into std::cerr if a cost or a constraint has been deleted.
     * In Release mode: No output.
     */
    void checkDeleteCostsAndConstraints();

protected:
    struct Constraints {
        Constraints();
        void clear();
        void updateNr();

        int nrEqConstr;
        int nrIneqConstr;
        std::vector<std::shared_ptr<Constraint>> spConstr;
        std::vector<std::shared_ptr<EqIneqConstraint>> spEqConstr;
        std::vector<std::shared_ptr<EqIneqConstraint>> spIneqConstr;
        std::vector<std::shared_ptr<ControlBoundConstraint>> spBoundConstr;
    };

protected:
    std::shared_ptr<PreviewSystem> ps_;
    std::unique_ptr<SolverInterface> sol_;
    std::vector<std::shared_ptr<CostFunction>> spCost_;
    Constraints constraints_;

    Eigen::MatrixXd Q_, Aineq_, Aeq_;
    Eigen::VectorXd c_, bineq_, beq_, lb_, ub_;

    std::chrono::duration<double> solveTime_, solveAndBuildTime_;
};

} // namespace pc
