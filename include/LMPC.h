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

/*! \brief The Linear Model Predictive Controller (LMPC).
 * This class gives all the needed composants for performing a model preview control.
 * It solves:\n
 * \f$X = \Phi x_{0} + \Psi U + \Xi\f$, where \f$U\f$ is the optimization vector.
 * \note \f$X = [x_0^T x_1^T ... x_N^T]^T\f$ and \f$U = [u_0^T u_1^T ... u_{N-1}^T]^T\f$
 * where \f$N\f$ is the dimension of the system (the number of steps).
 * \warning This class waits for a discretized system ! Continuous systems are not implemented.
 */
class COPRA_DLLAPI LMPC {
public:
    /*! \brief Default constructor.
     * Initialize problem variables to default and get the desired solver.
     * You need to call initializeController before using the MPCTypeFull.
     * \param sFlag The flag corresponding to the desired solver
     */
    LMPC(SolverFlag sFlag = SolverFlag::DEFAULT);
    /*! \brief Constructor.
     * Initialize problem variables w.r.t. the PreviewSystem and get the desired solver.
     * \param ps A preview system to make a copy from
     * \param sFlag The flag corresponding to the desired solver
     */
    LMPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag = SolverFlag::DEFAULT);
    /*! \brief Default virtual destructor. */
    virtual ~LMPC() = default;
    /*! \brief Default move-constructor. */
    LMPC(LMPC&&) = default;
    /*! \brief Default move-assign operator. */
    LMPC& operator=(LMPC&&) = default;

    /*! \brief Select a solver.
     * It will load a solver with default values.
     * \see useSolver if you want to give to the MPC a solver with specific parameters.
     * \param flag The solver to use \see pc::SolverFlag
     */
    void selectQPSolver(SolverFlag flag);
    /*! \brief Make the MPC uses a user-defined qp solver.
     * \param solver The user-defined solver. It must inheritate from the SolverInterface
     */
    void useSolver(std::unique_ptr<SolverInterface>&& solver);
    /*! \brief Initialize the controller with regard to the preview system.
     * This function needs to be called each time the system dimension changes.
     * \param ps The preview system
     */
    void initializeController(const std::shared_ptr<PreviewSystem>& ps);
    /*! \brief Solve the system.
     * Fill Phi, Psi, xi in PreviewSystem.
     * Fill A, b in Constraints.
     * \return True if a solution has been found.
     */
    bool solve();
    /*! \brief Print information on the QP solver status. */
    void inform() const noexcept;
    /*! \brief Return the time (in s) needed to solve the qp problem. */
    double solveTime() const noexcept;
    /*! \brief Return the time (in s) needed to build and solve the qp problem. */
    double solveAndBuildTime() const noexcept;
    /*! \brief Add a cost function to the system.
     * \param costFun A cost type \see TrajectoryCost \see TargetCost
     * \see ControlCost \see MixedTrajectoryCost \see MixedTargetCost
     */
    void addCost(const std::shared_ptr<CostFunction>& costFun);
    /*! \brief Add a constraint to the system.
     * \param constr A constraint type \see TrajectoryConstrain \see ControlConstraint
     * \see TrajectoryBoundConstraint \see ControlBoundConstraint
     */
    void addConstraint(const std::shared_ptr<Constraint>& constr);
    /*! \brief Clear all costs. */
    void clearCosts() noexcept;
    /*! \brief Clear all constraints. */
    void clearConstraints() noexcept;
    /*! \brief Remove specified cost. */
    void removeCost(const std::shared_ptr<CostFunction>& costFun);
    /*! \brief Remove specified constraint. */
    void removeConstraint(const std::shared_ptr<Constraint>& constr);

    // Getter Functions
    /*! \brief Return the control vector \f$U\f$. */
    const Eigen::VectorXd& control() const noexcept { return control_; }
    /*! \brief Return the trajectory vector \f$X\f$. */
    const Eigen::VectorXd& trajectory() const noexcept { return trajectory_; }
    /*! \brief Return the number of equality constraint. */
    inline int nrEqConstr() { return constraints_.nrEqConstr; }
    /*! \brief Return the number of inequality constraint. */
    inline int nrIneqConstr() { return constraints_.nrIneqConstr; }
    /*! \brief Return cost matrix. */
    inline const Eigen::MatrixXd& Q() { return Q_; }
    /*! \brief Return cost vector. */
    inline const Eigen::VectorXd& c() { return c_; }
    /*! \brief Return inequality constraint matrix. */
    inline const Eigen::MatrixXd& Aineq() { return Aineq_; }
    /*! \brief Return inequality constraint vector. */
    inline const Eigen::VectorXd& bineq() { return bineq_; }
    /*! \brief Return equality constraint matrix. */
    inline const Eigen::MatrixXd& Aeq() { return Aeq_; }
    /*! \brief Return equality constraint vector. */
    inline const Eigen::VectorXd& beq() { return beq_; }
    /*! \brief Return lower bound constraint vector. */
    inline const Eigen::VectorXd& lb() { return lb_; }
    /*! \brief Return upper bound constraint vector. */
    inline const Eigen::VectorXd& ub() { return ub_; }

private:
    /*! \brief Append the constraint into the QP \see Constraints. */
    void addConstraintByType(const std::shared_ptr<Constraint>& constr);

    /*! \brief Clear Aeq, beq, Aineq, bineq, ub, lb. */
    virtual void clearConstraintMatrices();
    /*! \brief Resize Q, c size from PreviewSystem. */
    virtual void updateQPMatrixSize();
    /*! \brief Resize Aeq, beq, Aineq, bineq, ub, lb size from constraints. */
    virtual void updateConstraintMatrixSize();
    /**! \brief Update the system and its constraints.
     * Fill A, b in Constraints.
     */
    virtual void updateSystem();
    /*! \brief Generate the QP-like format matrices. */
    virtual void makeQPForm();
    /*! \brief Update trajectory and control. */
    virtual void updateResults();
    /*! \brief Check for costs and constraints to delete.
     * Check if a cost or a constraint still exist.
     * In Debug mode: Output into std::cerr if a cost or a constraint has been deleted.
     * In Release mode: No output.
     */
    void checkDeleteCostsAndConstraints();

protected:
    /*! \brief Nested representation of Constraints. */
    struct Constraints {
        /*! \brief Default constructor. */
        Constraints();
        /*! \brief Clear all constraints. */
        void clear();
        /*! \brief Update the number of constraints. */
        void updateNr();

        int nrEqConstr; /*!< Number of equality constraint */
        int nrIneqConstr; /*!< Number of inequality constraint */
        std::vector<std::shared_ptr<Constraint>> spConstr; /*!< Vector of all constraints */
        std::vector<std::shared_ptr<EqIneqConstraint>> spEqConstr; /*!< Vector of all equality constraints */
        std::vector<std::shared_ptr<EqIneqConstraint>> spIneqConstr; /*!< Vector of all inequality constraints */
        std::vector<std::shared_ptr<ControlBoundConstraint>> spBoundConstr; /*!< Vector of all bound constraints */
    };

protected:
    std::shared_ptr<PreviewSystem> ps_; /*!< Preview System to work with */
    std::unique_ptr<SolverInterface> sol_; /*!< Underlying QP solver */
    std::vector<std::shared_ptr<CostFunction>> spCost_; /*!< Vector of all costs */
    Constraints constraints_; /*!< Set of all constraints */

    Eigen::MatrixXd Q_; /*!< Cost matrix */
    Eigen::VectorXd c_; /*!< Cost vector */
    Eigen::MatrixXd Aineq_; /*!< Inequality constraint matrix */
    Eigen::VectorXd bineq_; /*!< Inequality constraint vector */
    Eigen::MatrixXd Aeq_; /*!< Equality constraint matrix */
    Eigen::VectorXd beq_; /*!< Equality constraint vector */
    Eigen::VectorXd lb_; /*!< Lower bound constraint vector */
    Eigen::VectorXd ub_; /*!< Upper bound constraint vector */
    Eigen::VectorXd trajectory_; /*!< Trajectory result */
    Eigen::VectorXd control_; /*!< Control result */

    std::chrono::duration<double> solveTime_; /*!< Solving time */
    std::chrono::duration<double> solveAndBuildTime_; /*!< Solving + building time */
};

} // namespace copra
