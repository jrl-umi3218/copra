/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "debugUtils.h"
#include "typedefs.h"
#include <Eigen/Core>
#include <memory>

namespace copra {

/*! \brief Abstract base class that represents cost functions.
 * Any derived class of this one can be added to the MPC \see MPC::addCost.
 * A cost function is written \f$\frac{1}{2} x^T Q x + c^T x\f$.
 * Or is written by the extended version \f$\frac{1}{2} x^T \tilde{Q} x + \tilde{c}^T x\f$
 * with \f$\tilde{Q}=[RR^T+EQ^{\minus 1}E^T, E; E^T, Q]\f$ and \f$\tilde{c}=[r; f]\f$
 */
class COPRA_DLLAPI CostFunction {
public:
    /*! Constructor.
     * \param name Name of the constraint
     */
    CostFunction(std::string&& name);
    /*! \brief Default virtual desturctor. */
    virtual ~CostFunction() = default;

    /*! \brief Generate the full size matrices
     * This allows one to give a cost constant matrix/vector for one step
     * and another cost matrix/vector with the full horizon.
     * It will expand the first matrix/vector to fit the full horizon.
     */
    virtual void autoSpan();

    /*! \brief Initialize the cost function.
     * \param ps The preview system.
     * \throw Throw an std::domain_error if M, N or p is of bad dimension.
     */
    virtual void initializeCost(const PreviewSystem& ps);

    /*! \brief Compute the cost value for a particular control input
     * \param control The control vector \f$U\f$.
     * \throw Throw an error if cost is not yet initialized.
     */
    virtual double getCostValue(const Eigen::VectorXd& control);

    /*! Compute \f$Q\f$ and \f$c\f$.
     * \param ps The preview system.
     */
    virtual void update(const PreviewSystem& ps) = 0;

    /*! \brief Set the weights of the system.
     * Perform a move semantic if a fullsize vector is given as rvalue (this is faster).
     * \param weights of to set.
     * \throw Throw a std::domain_error if weights is badly dimension.
     */
    template <typename TVec, typename = std::enable_if_t<!is_all_arithmetic<TVec>::value>>
    void weights(TVec&& weights)
    {
        if (weights.rows() == weights_.rows()) {
            weights_ = std::forward<TVec>(weights);
        } else if (weights_.rows() % weights.rows() == 0) {
            auto size = weights_.rows() / weights.rows();
            for (auto i = 0; i < size; ++i) {
                weights_.segment(i * weights.rows(), weights.rows()) = weights;
            }
        } else {
            DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnDim("weights", weights, weights_.rows()));
        }
    }
    /*! \brief Set the weights of the system.
     * All variables are set to the same weight.
     * \param weight Weight to apply.
     */
    template <typename T, typename = std::enable_if_t<is_all_arithmetic<T>::value>>
    void weight(T weight)
    {
        weights_.setConstant(weight);
    }

    /*! \brief Return the name of the constraint. */
    const std::string& name() const noexcept { return name_; }
    /*! \brief Return the state-vector cost matrix. */
    const Eigen::MatrixXd& Q() const noexcept { return Q_; }
    /*! \brief Return the state-vector cost vector. */
    const Eigen::VectorXd& c() const noexcept { return c_; }
    /*! \brief Return the cross intial state-vector/state-vector cost matrix. */
    const Eigen::MatrixXd& E() const noexcept { return E_; }
    /*! \brief Return the state-vector cost vector for extended version. */
    const Eigen::VectorXd& f() const noexcept { return f_; }

protected:
    std::string name_; /*!< Name of the cost function */
    bool fullSizeEntry_; /*!< False if step-vector/step-matrix, True otherwise */
    Eigen::MatrixXd Q_; /*!< State-vector cost matrix */
    Eigen::VectorXd c_; /*!< State-vector cost vector */
    Eigen::MatrixXd E_; /*!< Initial State-vector cost matrix */
    Eigen::VectorXd f_; /*!< Initial State-vector cost vector */
    Eigen::VectorXd weights_; /*!< Cost weights */
};

/*! \brief Trajectory cost function class.
 * This cost function looks for a minimization around a trajectory.
 * Mathematically, it is \f$(MX-p)^TW_X(MX-p) \Leftrightarrow \sum_k w_k\|Mx_k-p\|^2\f$.
 */
class COPRA_DLLAPI TrajectoryCost final : public CostFunction {
public:
    /*! \brief Constructor of the trajectory cost function.
     * Create a cost function of type \f$(MX-p)^TW_X(MX-p)\f$.
     * Perform a move semantic if an rvalue is given.
     * \param M The matrix side of the cost function
     * \param p The vector side of the cost function
     */
    template <typename TMat, typename TVec, typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value>>
    TrajectoryCost(TMat&& M, TVec&& p)
        : CostFunction("TrajectoryCost")
        , M_(std::forward<TMat>(M))
        , p_(std::forward<TVec>(p))
    {
        weights_ = Eigen::VectorXd::Ones(p_.rows());
    }

    void autoSpan() override;
    void initializeCost(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;

private:
    Eigen::MatrixXd M_; /*!< State-vector constraint matrix */
    Eigen::VectorXd p_; /*!< Bias constraint vector */
};

/*! \brief Target cost function class.
 * This cost function looks for target a final point.
 * Mathematically, it is \f$(Mx_N-p)^Tw_x(Mx_N-p) \Leftrightarrow w_x\|Mx_N-p\|^2\f$.
 */
class COPRA_DLLAPI TargetCost final : public CostFunction {
public:
    /*! \brief Constructor of the target cost function.
     * Create a cost function of type \f$(Mx_N-p)^TW_X(Mx_N-p)\f$.
     * Perform a move semantic if an rvalue is given.
     * \param M The matrix side of the cost function
     * \param p The vector side of the cost function
     */
    template <typename TMat, typename TVec, typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value>>
    TargetCost(TMat&& M, TVec&& p)
        : CostFunction("TargetCost")
        , M_(std::forward<TMat>(M))
        , p_(std::forward<TVec>(p))
    {
        weights_ = Eigen::VectorXd::Ones(p_.rows());
    }

    void initializeCost(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;

private:
    Eigen::MatrixXd M_; /*!< State-vector constraint matrix */
    Eigen::VectorXd p_; /*!< Bias constraint vector */
};

/*! \brief Control cost function class.
 * This cost function looks for a minimization of the control.
 * Mathematically, it is \f$(NU-p)^TW_U(NU-p) \Leftrightarrow sum_k w_u\|Nu_k-p\|^2\f$.
 */
class COPRA_DLLAPI ControlCost final : public CostFunction {
public:
    /*! \brief Constructor of the control cost function.
     * Create a cost function of type \f$(NU-p)^TW_U(NU-p)\f$.
     * Perform a move semantic if an rvalue is given.
     * \param N The matrix side of the cost function
     * \param p The vector side of the cost function
     */
    template <typename TMat, typename TVec, typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value>>
    ControlCost(TMat&& N, TVec&& p)
        : CostFunction("ControlCost")
        , N_(std::forward<TMat>(N))
        , p_(std::forward<TVec>(p))
    {
        weights_ = Eigen::VectorXd::Ones(p_.rows());
    }

    void autoSpan() override;
    void update(const PreviewSystem& ps) override;
    void initializeCost(const PreviewSystem& ps) override;

private:
    Eigen::MatrixXd N_; /*!< Control-vector constraint matrix */
    Eigen::VectorXd p_; /*!< Bias constraint vector */
};

/*! \brief Mixed cost function class.
 * This cost function looks for a minimization of a linear combination of trajectory and control.
 * Mathematically, it is \f$(MX+NU-p)^TW_M(MX+NU-p) \Leftrightarrow sum_k w_m\|Mx_k+Nu_k-p\|^2\f$.
 */
class COPRA_DLLAPI MixedCost final : public CostFunction {
public:
    /*! \brief Constructor of the mixed cost function.
     * Create a cost function of type \f$(MX+NU-p)^TW_M(MX+NU-p)\f$.
     * Perform a move semantic if an rvalue is given.
     * \param M The matrix side of the trajectory cost function
     * \param N The matrix side of the control cost function
     * \param p The vector side of the cost function
     */
    template <typename TMat1, typename TMat2, typename TVec, typename = std::enable_if_t<!is_all_arithmetic<TMat1, TMat2, TVec>::value>>
    MixedCost(TMat1&& M, TMat2&& N, TVec&& p)
        : CostFunction("MixedCost")
        , M_(std::forward<TMat1>(M))
        , N_(std::forward<TMat2>(N))
        , p_(std::forward<TVec>(p))
    {
        weights_ = Eigen::VectorXd::Ones(p_.rows());
    }

    void autoSpan() override;
    void update(const PreviewSystem& ps) override;
    void initializeCost(const PreviewSystem& ps) override;

private:
    Eigen::MatrixXd M_; /*!< State-vector constraint matrix */
    Eigen::MatrixXd N_; /*!< Control-vector constraint matrix */
    Eigen::VectorXd p_; /*!< Bias constraint vector */
};

} // namespace copra
