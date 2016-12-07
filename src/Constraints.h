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

#include <Eigen/Core>
#include <string>

namespace mpc {

//forward declaration
struct PreviewSystem;

/**
 * Flags to identify the type of the constraints
 */
enum class ConstraintFlag {
    EqualityConstraint, /**< Equality constraint flag */
    InequalityConstraint, /**< Inequality constraint flag */
    BoundConstraint /**< Bound constraint tag */
};

/**
 * Abstract base class that represents Constraints.
 * Any derived class of this one can be added to the PreviewController @see MPCTypeFull::addConstraint.
 */
class Constraint {
public:
    /**
     * Constructor of a constraint.
     * @param name Name of the constraint
     */
    Constraint(const std::string& name);

    /**
     * Declare virtual desturctor
     */
    virtual ~Constraint() = default;

    /**
     * Initialization of the constraint.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void initializeConstraint(const PreviewSystem& ps) = 0;

    /**
     * Update the constraint.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void update(const PreviewSystem& ps) = 0;

    /**
     * Get the type of the constraint
     * @return The type of the constraint @see ConstraintFlag
     */
    virtual ConstraintFlag constraintType() = 0;

    /**
     * Function that return the name of the constraint.
     * @return The name of the constraint.
     */
    const std::string& name() const noexcept
    {
        return name_;
    }

    /**
     * Function that return the number of constraints.
     * @return The number of constraints.
     */
    int nrConstr() noexcept
    {
        return nrConstr_;
    }

protected:
    std::string name_;
    int nrConstr_;
};

/**
 * Abstract Class for Equality and Inequality constraints.
 * Even if Equality and Inequality constraints are different,
 * their matrices are written the same way.
 */
class EqIneqConstraint : public Constraint {
public:
    /**
     * Constructor of a constraint.
     * @param name Name of the constraint
     * @param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     */
    EqIneqConstraint(const std::string& name, bool isInequalityConstraint);

    /**
     * Get the 'A' matrix of the Equality/Inequality (\f$Ax\leq b\f$, \f$Ax = b\f$)
     * @return The 'A' matrix of the constraint
     */
    const Eigen::MatrixXd& A()
    {
        return A_;
    }

    /**
     * Get the 'b' vector of the Equality/Inequality (\f$Ax\leq b\f$, \f$Ax = b\f$)
     * @return The 'b' vector of the constraint
     */
    const Eigen::VectorXd& b()
    {
        return b_;
    }

protected:
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    bool isIneq_;
};

/**
 * Trajectory constraint.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints \f$Ex = f\f$
 * or an Inequality constraints \f$Ex\leq f\f$.
 */
class TrajectoryConstraint final : public EqIneqConstraint {
public:
    /**
     * Constructor of the trajectory constraint.
     * Create a constraint of type \f$Ex\leq f\f$ or \f$Ex = f\f$.
     * As \f$U\f$ is the optimization variable, \f$Ex\leq f\f$ or \f$Ex = f\f$
     * is transformed to be \f$AU\leq b\f$ or \f$AU = b\f$.
     * @warning \f$E\f$ and \f$f\f$ size must be the size of one iteration.
     * @param E The matrix side of the constraint
     * @param f The vector side of the constraint
     * @param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     */
    TrajectoryConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint = true);

    /**
     * Initialize the constraint by resizing its inner matrices and vectors
     * and setting the number of constraints.
     * @param ps A preview system.
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$A\f$ and \f$b\f$ from \f$E\f$, \f$f\f$ and the preview system.
     * @param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * @return @see ConstraintFlag::InequalityConstraint if constructor's isInequalityConstraint is true
     * @return @see ConstraintFlag::EqualityConstraint if constructor's isInequalityConstraint is false
     */
    ConstraintFlag constraintType() override;

private:
    Eigen::MatrixXd E_;
    Eigen::VectorXd f_;
};

/**
 * Control constraint.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints \f$Eu = f\f$
 * or an Inequality constraints \f$Eu\leq f\f$.
 * @todo Gives the possiblity to set directly \f$A\f$ and \f$b\f$. 
 */
class ControlConstraint final : public EqIneqConstraint {
public:
    /**
     * Constructor of the control constraint.
     * Create a constraint of type \f$Eu\leq f\f$ or \f$Eu = f\f$.
     * As \f$U\f$ is the optimization variable, \f$Eu\leq f\f$ or \f$Eu = f\f$
     * is transformed to be \f$AU\leq b\f$ or \f$AU = b\f$.
     * @warning \f$E\f$ and \f$f\f$ size must be the size of one iteration.
     * @param E The matrix side of the constraint
     * @param f The vector side of the constraint
     * @param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     */
    ControlConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint = true);

    /**
     * Initialize the constraint by resizing its inner matrices and vectors
     * and setting the number of constraints.
     * @param ps A preview system.
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$A\f$ and \f$b\f$ from \f$E\f$, \f$f\f$ and the preview system.
     * @param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * @return @see ConstraintFlag::InequalityConstraint if constructor's isInequalityConstraint is true
     * @return @see ConstraintFlag::EqualityConstraint if constructor's isInequalityConstraint is false
     */
    ConstraintFlag constraintType() override;

private:
    Eigen::MatrixXd E_;
    Eigen::VectorXd f_;
    bool isIneq_;
};

/**
 * Trajectory Bound constraint.
 * Even if it is a bound constraint, the optimization vector is \f$U\f$ 
 * so this constraint has to be transformed to an Inequality constraint.
 * @warning This constraint is defined in the QP as an Inequality constraint. 
 */
class TrajectoryBoundConstraint final : public EqIneqConstraint {
public:
    /**
     * Constructor of the trajectory Bound constraint.
     * Create a constraint of type \f$\underbar{x}\leq x\leq\bar{x}\f$.
     * As \f$U\f$ is the optimization variable,
     * it is transformed to be \f$AU\leq b\f$.
     * @warning \f$\underbar{x}\f$ and \f$\bar{x}\f$ size must be the size of one iteration.
     * @param lower The lower bound \f$\underbar{x}\f$ of the constraint
     * @param upper The upper bound \f$\bar{x}\f$ of the constraint
     */
    TrajectoryBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper);

    /**
     * Initialize the constraint by resizing its inner matrices and vectors
     * and setting the number of constraints.
     * @param ps A preview system.
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$A\f$ and \f$b\f$ from \f$E\f$, \f$f\f$ and the preview system.
     * @param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * @return @see ConstraintFlag::InequalityConstraint
     */
    ConstraintFlag constraintType() override;

private:
    Eigen::VectorXd lower_, upper_;
};

/**
 * Control Bound constraint.
 * It bounds the optimization \f$\underbar{u}\leq u \leq\bar{u}\f$
 * @todo Gives the possiblity to set directly \f$\underbar{U}\f$ and \f$\bar{U}\f$. 
 */
class ControlBoundConstraint final : public Constraint {
public:
    /**
     * Constructor of the trajectory Bound constraint.
     * Create a constraint of type \f$\underbar{u}\leq u\leq\bar{u}\f$.
     * As \f$U\f$ is the optimization variable,
     * it is transformed to be \f$\underbar{U}\leq U\leq\bar{U}\f$.
     * @warning \f$\underbar{u}\f$ and \f$\bar{u}\f$ size must be the size of one iteration.
     * @param lower The lower bound \f$\underbar{u}\f$ of the constraint
     * @param upper The upper bound \f$\bar{u}\f$ of the constraint
     */
    ControlBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper);

    /**
     * Initialize the constraint by resizing its inner matrices and vectors
     * and setting the number of constraints.
     * @param ps A preview system.
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$\underbar{U}\f$ and \f$\bar{U}\f$ from \f$\underbar{u}\f$, \f$\bar{u}\f$ and the preview system.
     * @param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * @return @see ConstraintFlag::BoundConstraint
     */
    ConstraintFlag constraintType() override;

    const Eigen::VectorXd& lower()
    {
        return lb_;
    }

    const Eigen::VectorXd& upper()
    {
        return ub_;
    }

private:
    Eigen::VectorXd lower_, upper_, lb_, ub_;
};

} // namespace mpc