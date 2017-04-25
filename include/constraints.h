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
#include <string>
#include <vector>

// eigen
#include <Eigen/Core>

// mpc
#include "config.hh"
#include "debugUtils.h"
#include "typedefs.h"

namespace mpc {

//forward declaration
struct PreviewSystem;

/**
 * Flags to identify the type of the constraints
 */
enum class ConstraintFlag {
    Constraint, /**< Any constraints */
    EqualityConstraint, /**< Equality constraint flag */
    InequalityConstraint, /**< Inequality constraint flag */
    BoundConstraint /**< Bound constraint tag */
};

/**
 * \brief Abstract base class that represents Constraints.
 * Any derived class of this one can be added to the PreviewController \see MPCTypeFull::addConstraint.
 */
class MPC_DLLAPI Constraint {
public:
    /**
     * Constructor of a constraint.
     * \param name Name of the constraint
     */
    Constraint(std::string&& name);

    /**
     * Declare virtual desturctor
     */
    virtual ~Constraint() = default;

    /**
     * \brief Generate the full size matrices
     * This allows one to give a constraint constant matrix/vector for one step
     * and another constraint matrix/vector with the full horizon.
     * It will expand the first matrix/vector to fit the full horizon.
     */
    virtual void autoSpan() = 0;

    /**
     * Initialization of the constraint.
     * \param ps An unique pointer to a PreviewSystem.
     */
    virtual void initializeConstraint(const PreviewSystem& ps) = 0;

    /**
     * Update the constraint.
     * \param ps An unique pointer to a PreviewSystem.
     */
    virtual void update(const PreviewSystem& ps) = 0;

    /**
     * Get the type of the constraint
     * \return The type of the constraint \see ConstraintFlag
     */
    virtual ConstraintFlag constraintType() const noexcept = 0;

    /**
     * Function that return the name of the constraint.
     * \return The name of the constraint.
     */
    const std::string& name() const noexcept
    {
        return name_;
    }

    /**
     * Function that return the number of constraints.
     * \return The number of constraints.
     */
    int nrConstr() noexcept
    {
        return nrConstr_;
    }

protected:
    std::string name_;
    int nrConstr_;
    bool fullSizeEntry_;
    bool hasBeenInitialized_;
};

/**
 * \brief Abstract Class for Equality and Inequality constraints.
 * Even if Equality and Inequality constraints are different, their matrices are written the same way.
 */
class MPC_DLLAPI EqIneqConstraint : public Constraint {
public:
    /**
     * Constructor of a constraint.
     * \param name Name of the constraint
     * \param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     */
    EqIneqConstraint(const std::string& name, bool isInequalityConstraint);

    /**
     * Get the 'A' matrix of the Equality/Inequality (\f$Ax\leq b\f$, \f$Ax = b\f$)
     * \return The 'A' matrix of the constraint
     */
    const Eigen::MatrixXd& A() const noexcept
    {
        return A_;
    }

    /**
     * Get the 'b' vector of the Equality/Inequality (\f$Ax\leq b\f$, \f$Ax = b\f$)
     * \return The 'b' vector of the constraint
     */
    const Eigen::VectorXd& b() const noexcept
    {
        return b_;
    }

protected:
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    bool isIneq_;
};

/**
 * \brief Trajectory constraint class.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints (\f$Ex = f\f$ or \f$EX = f\f$)
 * or an Inequality constraints (\f$Ex\leq f\f$ or \f$EX\leq f\f$) with \f$X=[x_1^T ... x_{nrStep}^T]^T\f$.
 */
class MPC_DLLAPI TrajectoryConstraint final : public EqIneqConstraint {
public:
    /**
     * \brief Constructor of the trajectory constraint.
     * Create a constraint of type \f$Ex\leq f\f$ or \f$Ex = f\f$ or \f$EX\leq f\f$ or \f$EX = f\f$ with \f$X=[x_1^T ... x_N^T]^T\f$.\n
     * As \f$U\f$ is the optimization variable, \f$Ex\leq f\f$ or \f$Ex = f\f$
     * is transformed to be \f$AU\leq b\f$ or \f$AU = b\f$.
     * Perform a move semantic if an rvalue is given (this is faster).
     * \param E The matrix side of the constraint
     * \param f The vector side of the constraint
     * \param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false)
     * \throw Throw an std::domain_error if E and f have not the same number of rows
     */
    template <typename TMat, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value> >
    TrajectoryConstraint(TMat&& E, TVec&& f, bool isInequalityConstraint = true)
        : EqIneqConstraint("Trajectory", isInequalityConstraint)
        , E_(std::forward<TMat>(E))
        , f_(std::forward<TVec>(f))
    {
    }

    /**
     * \brief Generate the full size matrices
     * If you have create the constraint with matrix \f$E_k\f$ and vector \f$h\f$
     * or \f$E\f$ and \f$h_k\f$ you need to call this function to redimension the matrxix/vector.
     */
    void autoSpan() override;

    /**
     * \brief Initialize the constraint.
     * This is done by resizing its inner matrices and vectors and setting the number of constraints.
     * \param ps A preview system.
     * \throw Throw an std::domain_error if E or f is not of the dimension of the preview system.
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$A\f$ and \f$b\f$ from \f$E\f$, \f$f\f$ and the preview system.
     * \param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * \return \see ConstraintFlag::InequalityConstraint if constructor's isInequalityConstraint is true
     * \return \see ConstraintFlag::EqualityConstraint if constructor's isInequalityConstraint is false
     */
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::MatrixXd E_;
    Eigen::VectorXd f_;
};

/**
 * \brief Control constraint class.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints (\f$Eu = f\f$ or \f$EU = f\f$)
 * or an Inequality constraints (\f$Eu\leq f\f$ or \f$EU\leq f\f$) with \f$U=[u_0^T ... u_{N-1}^T]^T\f$
 */
class MPC_DLLAPI ControlConstraint final : public EqIneqConstraint {
public:
    /**
     * \brief Constructor of the control constraint.
     * Create a constraint of type \f$Gu\leq f\f$ or \f$Gu = f\f$ or \f$GU = f\f$ or \f$GU\leq f\f$ with \f$U=[u_0^T ... u_{N-1}^T]^T\f$.\n
     * As \f$U\f$ is the optimization variable, \f$Gu\leq f\f$ or \f$Gu = f\f$
     * is transformed to be \f$AU\leq b\f$ or \f$AU = b\f$.
     * Perform a move semantic if an rvalue is given (this is faster).
     * \param G The matrix side of the constraint
     * \param f The vector side of the constraint
     * \param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     * \throw Throw an std::domain_error if G and f have not the same number of rows
     */
    template <typename TMat, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value> >
    ControlConstraint(TMat&& G, TVec&& f, bool isInequalityConstraint = true)
        : EqIneqConstraint("Control", isInequalityConstraint)
        , G_(std::forward<TMat>(G))
        , f_(std::forward<TVec>(f))
    {
    }

    /**
     * \brief Generate the full size matrices
     * If you have create the constraint with matrix \f$G_k\f$ and vector \f$h\f$
     * or \f$G\f$ and \f$h_k\f$ you need to call this function to redimension the matrxix/vector.
     */
    void autoSpan() override;

    /**
     * \brief Initialize the constraint.
     * This is done by resizing its inner matrices and vectors and setting the number of constraints.
     * \param ps A preview system.
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$A\f$ and \f$b\f$ from \f$E\f$, \f$f\f$ and the preview system.
     * \param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * \return \see ConstraintFlag::InequalityConstraint if constructor's isInequalityConstraint is true
     * \return \see ConstraintFlag::EqualityConstraint if constructor's isInequalityConstraint is false
     */
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::MatrixXd G_;
    Eigen::VectorXd f_;
};

/**
 * \brief Mixed constraint class.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints (\f$Ex + Gu = f\f$ or \f$EX + GU = f\f$)\n
 * or an Inequality constraints (\f$Ex + Gu\leq f\f$ or \f$EX + GU\leq f\f$)\n
 * with \f$X=[x_1^T ... x_N^T]^T\f$ and \f$U=[u_0^T ... u_{N-1}^T]^T\f$
 */
class MPC_DLLAPI MixedConstraint final : public EqIneqConstraint {
public:
    /**
     * \brief Constructor of the control constraint.
     * Create a constraint of type \f$Ex + Gu\leq f\f$ or \f$Ex + Gu = f\f$ or \f$EX + GU = f\f$ or \f$EX + GU\leq f\f$\n
     * with \f$X=[x_1^T ... x_N^T]^T\f$ and \f$U=[u_0^T ... u_{N-1}^T]^T\f$.\n
     * As \f$U\f$ is the optimization variable, \f$Ex + Gu\leq f\f$ or \f$Ex + Gu = f\f$
     * is transformed to be \f$AU\leq b\f$ or \f$AU = b\f$.
     * Perform a move semantic if an rvalue is given (this is faster).
     * \note Please use \see ControlConstraint and \see TrajectoryConstraint for non-mixed constraint (they are sligthly faster)
     * \param E The matrix applied to the trajectory part of the constraint
     * \param G The matrix applied to the control part of the constraint
     * \param f The vector side of the constraint
     * \param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     * \throw Throw an std::domain_error if E, G and f have not the same number of rows
     */
    template <typename TMat1, typename TMat2, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat1, TMat2, TVec>::value> >
    MixedConstraint(TMat1&& E, TMat2&& G, TVec&& f, bool isInequalityConstraint = true)
        : EqIneqConstraint("Control", isInequalityConstraint)
        , E_(std::forward<TMat1>(E))
        , G_(std::forward<TMat2>(G))
        , f_(std::forward<TVec>(f))
    {
    }

    /**
     * \brief Generate the full size matrices
     * If you have create the constraint with matrix \f$E_k\f$ \f$G_k\f$ and vector \f$h\f$
     * or \f$E_k\f$, \f$G\f$ and \f$h_k\f$, etc... you need to call this function to redimension the matrxix/vector.
     */
    void autoSpan() override;

    /**
     * \brief Initialize the constraint.
     * This is done by resizing its inner matrices and vectors and setting the number of constraints.
     * \param ps A preview system.
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$A\f$ and \f$b\f$ from \f$E\f$, \f$G\f$, \f$f\f$ and the preview system.
     * \param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * \return \see ConstraintFlag::InequalityConstraint if constructor's isInequalityConstraint is true
     * \return \see ConstraintFlag::EqualityConstraint if constructor's isInequalityConstraint is false
     */
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::MatrixXd E_, G_;
    Eigen::VectorXd f_;
};

/**
 * \brief Trajectory Bound constraint.
 * Even if it is a bound constraint, the optimization vector is \f$U\f$ 
 * so this constraint has to be transformed to an Inequality constraint.
 * \warning This constraint is defined in the QP as an Inequality constraint. 
 * It might be faster to transform yourself this constraint into an inequality constraint.
 */
class MPC_DLLAPI TrajectoryBoundConstraint final : public EqIneqConstraint {
public:
    /**
     * \brief Constructor of the trajectory Bound constraint.
     * Create a constraint of type \f$\underline{x}\leq x\leq\overline{x}\f$ or \f$\underline{X}\leq X\leq\overline{X}\f$ 
     * with \f$X=[x_1^T ... x_N^T]^T\f$.\n
     * As \f$U\f$ is the optimization variable,
     * \f$\underline{x}\leq x\leq\overline{x}\f$ is transformed to be \f$AU\leq b\f$.
     * Perform a move semantic if an rvalue is given (this is faster).
     * \param lower The lower bound \f$\underline{x}\f$ of the constraint
     * \param upper The upper bound \f$\overline{x}\f$ of the constraint
     * \throw Throw an std::domain_error if lower and upper are not of the same dimension
     */
    template <typename TVec1, typename TVec2,
        typename = std::enable_if_t<!is_all_arithmetic<TVec1, TVec2>::value> >
    TrajectoryBoundConstraint(TVec1&& lower, TVec2&& upper)
        : EqIneqConstraint("Trajectory bound", true)
        , lower_(std::forward<TVec1>(lower))
        , upper_(std::forward<TVec2>(upper))
        , lowerLines_()
        , upperLines_()
    {
        for (auto line = 0; line < lower_.rows(); ++line) {
            if (lower_(line) != -std::numeric_limits<double>::infinity())
                lowerLines_.push_back(line);
        }
        for (auto line = 0; line < upper_.rows(); ++line) {
            if (upper_(line) != std::numeric_limits<double>::infinity())
                upperLines_.push_back(line);
        }
    }

    /**
     * \brief Generate the full size matrices
     * If you have create the constraint with vector \f$\underline{x}\f$ and vector \f$\overline{X}\f$
     * or \f$\underline{X}\f$ and \f$\overline{x}\f$ you need to call this function to redimension the vectors.
     */
    void autoSpan() override;

    /**
     * \brief Initialize the constraint.
     * This is done by resizing its inner matrices and vectors and setting the number of constraints.
     * \param ps A preview system.
     * \throw Throw an std::domain_error if lower or upper is badly dimension
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$A\f$ and \f$b\f$ from \f$E\f$, \f$f\f$ and the preview system.
     * \param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * \return \see ConstraintFlag::InequalityConstraint
     */
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::VectorXd lower_, upper_;
    std::vector<int> lowerLines_, upperLines_;
};

/**
 * \brief Control Bound constraint.
 * It bounds the optimization \f$\underline{u}\leq u \leq\overline{u}\f$
 */
class MPC_DLLAPI ControlBoundConstraint final : public Constraint {
public:
    /**
     * \brief Constructor of the trajectory Bound constraint.
     * Create a constraint of type \f$\underline{u}\leq u\leq\overline{u}\f$ or \f$\underline{U}\leq U\leq\overline{U}\f$ 
     * with \f$U=[u_0^T ... u_{N-1}^T]^T\f$.\n
     * As \f$U\f$ is the optimization variable,
     * \f$\underline{u}\leq u\leq\overline{u}\f$ is transformed to be \f$\underline{U}\leq U\leq\overline{U}\f$.
     * Perform a move semantic if an rvalue is given (this is faster).
     * \param lower The lower bound \f$\underline{u}\f$ of the constraint
     * \param upper The upper bound \f$\overline{u}\f$ of the constraint
     * \throw Throw an std::domain_error if lower and upper are not of the same dimension
     */
    template <typename TVec1, typename TVec2,
        typename = std::enable_if_t<!is_all_arithmetic<TVec1, TVec2>::value> >
    ControlBoundConstraint(TVec1&& lower, TVec2&& upper)
        : Constraint("Control bound constraint")
        , lower_(std::forward<TVec1>(lower))
        , upper_(std::forward<TVec2>(upper))
        , lb_()
        , ub_()
    {
    }

    /**
     * \brief Generate the full size matrices
     * If you have create the constraint with vector \f$\underline{u}\f$ and vector \f$\overline{U}\f$
     * or \f$\underline{U}\f$ and \f$\overline{u}\f$ you need to call this function to redimension the vectors.
     */
    void autoSpan() override;

    /**
     * \brief Initialize the constraint.
     * This is by resizing its inner matrices and vectors and setting the number of constraints.
     * \param ps A preview system.
     * \throw Throw an std::domain_error if lower or upper is badly dimension
     */
    void initializeConstraint(const PreviewSystem& ps) override;

    /**
     * Compute \f$\underline{U}\f$ and \f$\overline{U}\f$ from \f$\underline{u}\f$, \f$\overline{u}\f$ and the preview system.
     * \param ps A preview system.
     */
    void update(const PreviewSystem& ps) override;

    /**
     * Get the type of the constraint
     * \return \see ConstraintFlag::BoundConstraint
     */
    ConstraintFlag constraintType() const noexcept override;

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
