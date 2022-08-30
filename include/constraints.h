/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "debugUtils.h"
#include "typedefs.h"
#include <Eigen/Core>
#include <string>
#include <vector>

namespace copra {

//forward declaration
struct PreviewSystem;

#if defined(__GNUC__)
#pragma GCC diagnostic push
// Work around GCC (< 6) bug see: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=43407
// The error still exist for GCC 7.4.0
#pragma GCC diagnostic ignored "-Wattributes"
#endif

/*! \brief Flags to identify the type of the constraints. */
enum class COPRA_DLLAPI ConstraintFlag {
    Constraint, /*!< Any constraints */
    EqualityConstraint, /*!< Equality constraint flag */
    InequalityConstraint, /*!< Inequality constraint flag */
    BoundConstraint /*!< Bound constraint tag */
};

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

/*! \brief Abstract base class that represents constraints.
 * Any derived class of this one can be added to the MPC \see MPC::addConstraint.
 */
class COPRA_DLLAPI Constraint {
public:
    /*! \brief Default constructor.
     * \param name Name of the constraint
     */
    Constraint(std::string&& name);
    /*! \brief Virtual destructor. */
    virtual ~Constraint() = default;

    /*! \brief Generate the full size matrices.
     * This allows one to give a constraint constant matrix/vector for one step
     * and another constraint matrix/vector with the full horizon.
     * It will expand the first matrix/vector to fit the full horizon.
     */
    virtual void autoSpan() = 0;
    /*! \brief Initialization of the constraint from the preview system.
     * \param ps The PreviewSystem
     * \throw Can throw an std::domain_error the constraint matrix size differ from the preview system.
     */
    virtual void initializeConstraint(const PreviewSystem& ps) = 0;
    /*! \brief Update the constraint from the preview system.
     * Compute \f$Y\f$, \f$A\f$, \f$z\f$ and \f$b\f$ from the constraint and the preview system.
     * \param ps The PreviewSystem
     */
    virtual void update(const PreviewSystem& ps) = 0;
    /*! \brief Check if the constraint is satisfied for a particular control vector.
     * \param control The control vector
     * \param tolerance The tolerance
     */
    virtual bool isSatisfied(const Eigen::VectorXd& control, const double tolerance=0) = 0;
    /*! \brief Return the type of the constraint. \see ConstraintFlag */
    virtual ConstraintFlag constraintType() const noexcept = 0;
    /*! \brief Return the name of the constraint. */
    const std::string& name() const noexcept { return name_; }
    /*! \brief Return the number of constraints. */
    int nrConstr() noexcept { return nrConstr_; }

protected:
    std::string name_; /*!< Name of the constraint */
    int nrConstr_; /*!< Number of constraints */
    bool fullSizeEntry_; /*!< False if step-vector/step-matrix, True otherwise. */
    bool hasBeenInitialized_; /*!< True if constraint already initialized, False otherwise. */
};

/*! \brief Abstract Class for Equality and Inequality constraints.
 * Even if Equality and Inequality constraints are different, their matrices are written the same way.
 */
class COPRA_DLLAPI EqIneqConstraint : public Constraint {
public:
    /*! \brief Constructor.
     * \param name Name of the constraint
     * \param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     */
    EqIneqConstraint(const std::string& name, bool isInequalityConstraint);

    /*! \brief Return the 'A' matrix of the Equality/Inequality (\f$Ax =|\leq b\f$ or \f$[Y, A]x =|\leq [z^T b^T]^T\f$). */
    const Eigen::MatrixXd& A() const noexcept { return A_; }
    /*! \brief Return the 'b' vector of the Equality/Inequality (\f$Ax = b\f$ or \f$[Y, A]x = [z^T b^T]^T\f$). */
    const Eigen::VectorXd& b() const noexcept { return b_; }
    /*! \brief Return the 'Y' matrix of the Equality/Inequality (\f$[Y, A]x = [z^T b^T]^T\f$). */
    const Eigen::MatrixXd& Y() const noexcept { return Y_; }
    /*! \brief Return the 'z' vector of the Equality/Inequality (\f$[Y, A]x = [z^T b^T]^T\f$). */
    const Eigen::VectorXd& z() const noexcept { return z_; }

protected:
    Eigen::MatrixXd A_; /*!< Constraint matrix of the state-vcector x */
    Eigen::VectorXd b_; /*!< Constraint vector of the state-vcector x */
    Eigen::MatrixXd Y_; /*!< Constraint matrix of the initial state */
    Eigen::VectorXd z_; /*!< Constraint vector of the initial state */
    bool isIneq_; /*!< True if inequality, False if equality */
};

/*! \brief Trajectory constraint class.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints (\f$Ex = f\f$ or \f$EX = f\f$)
 * or an Inequality constraints (\f$Ex\leq f\f$ or \f$EX\leq f\f$) with \f$X=[x_1^T ... x_{nrStep}^T]^T\f$.
 */
class COPRA_DLLAPI TrajectoryConstraint final : public EqIneqConstraint {
public:
    /*! \brief Constructor of the trajectory constraint.
     * Create a constraint of type \f$Ex\leq f\f$ or \f$Ex = f\f$ or \f$EX\leq f\f$ or \f$EX = f\f$ with \f$X=[x_1^T ... x_N^T]^T\f$.\n
     * As \f$U\f$ is the optimization variable, \f$Ex\leq f\f$ or \f$Ex = f\f$
     * is transformed to be \f$AU\leq b\f$ or \f$AU = b\f$.
     * Perform a move semantic if an rvalue is given (this is faster).
     * \param E The matrix side of the constraint
     * \param f The vector side of the constraint
     * \param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false)
     */
    template <typename TMat, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value>>
    TrajectoryConstraint(TMat&& E, TVec&& f, bool isInequalityConstraint = true)
        : EqIneqConstraint("Trajectory", isInequalityConstraint)
        , E_(std::forward<TMat>(E))
        , f_(std::forward<TVec>(f))
    {
    }

    void autoSpan() override;
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    bool isSatisfied(const Eigen::VectorXd& control, const double tolerance=0) override;
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::MatrixXd E_; /*!< State-vector constraint matrix */
    Eigen::VectorXd f_; /*!< Bias constraint bias vector */
};

/*! \brief Control constraint class.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints (\f$Eu = f\f$ or \f$EU = f\f$)
 * or an Inequality constraints (\f$Eu\leq f\f$ or \f$EU\leq f\f$) with \f$U=[u_0^T ... u_{N-1}^T]^T\f$
 */
class COPRA_DLLAPI ControlConstraint final : public EqIneqConstraint {
public:
    /*! \brief Constructor of the control constraint.
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
        typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value>>
    ControlConstraint(TMat&& G, TVec&& f, bool isInequalityConstraint = true)
        : EqIneqConstraint("Control", isInequalityConstraint)
        , G_(std::forward<TMat>(G))
        , f_(std::forward<TVec>(f))
    {
    }

    void autoSpan() override;
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    bool isSatisfied(const Eigen::VectorXd& control, const double tolerance=0) override;
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::MatrixXd G_; /*!< Control-vector constraint matrix */
    Eigen::VectorXd f_; /*!< Bias constraint vector */
};

/*! \brief Mixed constraint class.
 * Depending on the parameter 'isInequalityConstraint' during the construction
 * it can be an Equality constraints (\f$Ex + Gu = f\f$ or \f$EX + GU = f\f$)\n
 * or an Inequality constraints (\f$Ex + Gu\leq f\f$ or \f$EX + GU\leq f\f$)\n
 * with \f$X=[x_1^T ... x_N^T]^T\f$ and \f$U=[u_0^T ... u_{N-1}^T]^T\f$
 */
class COPRA_DLLAPI MixedConstraint final : public EqIneqConstraint {
public:
    /*! \brief Constructor of the control constraint.
     * Create a constraint of type \f$Ex + Gu\leq f\f$ or \f$Ex + Gu = f\f$ or \f$EX + GU = f\f$ or \f$EX + GU\leq f\f$\n
     * with \f$X=[x_1^T ... x_N^T]^T\f$ and \f$U=[u_0^T ... u_{N-1}^T]^T\f$.\n
     * As \f$U\f$ is the optimization variable, \f$Ex + Gu\leq f\f$ or \f$Ex + Gu = f\f$
     * is transformed to be \f$AU\leq b\f$ or \f$AU = b\f$.
     * Perform a move semantic if an rvalue is given (this is faster).
     * \note Please use \see ControlConstraint and \see TrajectoryConstraint for non-mixed constraint (they are slightly faster)
     * \param E The matrix applied to the trajectory part of the constraint
     * \param G The matrix applied to the control part of the constraint
     * \param f The vector side of the constraint
     * \param isInequalityConstraint Whether the constraint is an Inequality (true) or an Equality (false).
     * \throw Throw an std::domain_error if E, G and f have not the same number of rows
     */
    template <typename TMat1, typename TMat2, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat1, TMat2, TVec>::value>>
    MixedConstraint(TMat1&& E, TMat2&& G, TVec&& f, bool isInequalityConstraint = true)
        : EqIneqConstraint("Control", isInequalityConstraint)
        , E_(std::forward<TMat1>(E))
        , G_(std::forward<TMat2>(G))
        , f_(std::forward<TVec>(f))
    {
    }

    void autoSpan() override;
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    bool isSatisfied(const Eigen::VectorXd& control, const double tolerance=0) override;
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::MatrixXd E_; /*!< State-vector constraint matrix */
    Eigen::MatrixXd G_; /*!< Control-vector constraint matrix */
    Eigen::VectorXd f_; /*!< Bias constraint vector */
};

/*! \brief Trajectory Bound constraint.
 * Even if it is a bound constraint, the optimization vector is \f$U\f$ 
 * so this constraint has to be transformed to an Inequality constraint.
 * \warning This constraint is defined in the QP as an Inequality constraint. 
 * It might be faster to transform yourself this constraint into an inequality constraint.
 */
class COPRA_DLLAPI TrajectoryBoundConstraint final : public EqIneqConstraint {
public:
    /*! \brief Constructor of the trajectory bound constraint.
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
        typename = std::enable_if_t<!is_all_arithmetic<TVec1, TVec2>::value>>
    TrajectoryBoundConstraint(TVec1&& lower, TVec2&& upper)
        : EqIneqConstraint("Trajectory bound", true)
        , lower_(std::forward<TVec1>(lower))
        , upper_(std::forward<TVec2>(upper))
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

    void autoSpan() override;
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    bool isSatisfied(const Eigen::VectorXd& control, const double tolerance=0) override;
    ConstraintFlag constraintType() const noexcept override;

private:
    Eigen::VectorXd lower_; /*!< State-vector lower bound constraint vector */
    Eigen::VectorXd upper_; /*!< State-vector upper bound constraint vector */
    std::vector<int> lowerLines_; /*!< Set of lower constraints different from \f$-\infty\f$ */
    std::vector<int> upperLines_; /*!< Set of upper constraints different from \f$\infty\f$ */
};

/*! \brief Control Bound constraint.
 * It bounds the optimization \f$\underline{u}\leq u \leq\overline{u}\f$
 */
class COPRA_DLLAPI ControlBoundConstraint final : public Constraint {
public:
    /*! \brief Constructor of the control bound constraint.
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
        typename = std::enable_if_t<!is_all_arithmetic<TVec1, TVec2>::value>>
    ControlBoundConstraint(TVec1&& lower, TVec2&& upper)
        : Constraint("Control bound constraint")
        , lower_(std::forward<TVec1>(lower))
        , upper_(std::forward<TVec2>(upper))
    {
    }

    void autoSpan() override;
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    bool isSatisfied(const Eigen::VectorXd& control, const double tolerance=0) override;
    ConstraintFlag constraintType() const noexcept override;
    /*! \brief Return the lower bound of the constraint */
    const Eigen::VectorXd& lower() { return lb_; }
    /*! \brief Return the upper bound of the constraint */
    const Eigen::VectorXd& upper() { return ub_; }

private:
    Eigen::VectorXd lower_; /*!< Control-vector lower bound */
    Eigen::VectorXd upper_; /*!< Control-vector upper bound */
    Eigen::VectorXd lb_; /*!< Full control-vector lower bound */
    Eigen::VectorXd ub_; /*!< Full control-vector upper bound */
};

} // namespace copra
