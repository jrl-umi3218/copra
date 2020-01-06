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

/**
 * \brief Abstract base class that represents cost functions.
 * Any derived class of this one can be added to the MPC \see MPC::addCost.
 */
class COPRA_DLLAPI CostFunction {
public:
    /**
     * Constructor of a cost function.
     * \param name Name of the constraint
     */
    CostFunction(std::string&& name);

    /**
     * Declare virtual desturctor
     */
    virtual ~CostFunction() = default;

    /**
     * \brief Generate the full size matrices
     * This allows one to give a cost constant matrix/vector for one step
     * and another cost matrix/vector with the full horizon.
     * It will expand the first matrix/vector to fit the full horizon.
     */
    virtual void autoSpan();

    /**
     * \brief Initialize the cost function.
     * \param ps The preview system.
     * \throw Throw an std::domain_error if M, N or p is of bad dimension.
     */
    virtual void initializeCost(const PreviewSystem& ps);

    /**
     * Compute \f$Q\f$ and \f$c\f$.
     * \param ps The preview system.
     */
    virtual void update(const PreviewSystem& ps) = 0;

    /**
     * Set the weights of the system.
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
            for (auto i = 0; i < size; ++i)
                weights_.segment(i * weights.rows(), weights.rows()) = weights;
        } else {
            DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnDim("weights", weights, weights_.rows()));
        }
    }

    /**
     * Set the weights of the system. All variables are set to the same weight.
     * \param weight Weight to apply.
     */
    template <typename T, typename = std::enable_if_t<is_all_arithmetic<T>::value>>
    void weight(T weight)
    {
        weights_.setConstant(weight);
    }

    /**
     * Function that return the name of the cost function.
     * \return The name of the constraint.
     */
    const std::string& name() const noexcept { return name_; }

    /**
     * \brief Function that return the Q matrix of the cost function
     * A cost function is written \f$ \frac{1}{2} x^T Q x + c^T x\f$.
     * \return The Q matrix
     */
    const Eigen::MatrixXd& Q() const noexcept { return Q_; }

    /**
     * \brief Function that return the c vector of the cost function
     * A cost function is written \f$ \frac{1}{2} x^T Q x + c^T x\f$.
     * \return The c vector matrix
     */
    const Eigen::VectorXd& c() const noexcept { return c_; }

protected:
    std::string name_;
    bool fullSizeEntry_;
    Eigen::MatrixXd Q_;
    Eigen::VectorXd c_;
    Eigen::VectorXd weights_;
};

/**
 * \brief Trajectory cost function class.
 * This cost function looks for a minimization around a trajectory.
 * Mathematically, it is \f$(MX-p)^TW_X(MX-p) \Leftrightarrow \sum_k w_k\|Mx_k-p\|^2\f$.
 */
class COPRA_DLLAPI TrajectoryCost final : public CostFunction {
public:
    /**
     * \brief Constructor of the trajectory cost function.
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
    Eigen::MatrixXd M_;
    Eigen::VectorXd p_;
};

/**
 * \brief Target cost function class.
 * This cost function looks for target a final point.
 * Mathematically, it is \f$(Mx_N-p)^Tw_x(Mx_N-p) \Leftrightarrow w_x\|Mx_N-p\|^2\f$.
 */
class COPRA_DLLAPI TargetCost final : public CostFunction {
public:
    /**
     * \brief Constructor of the target cost function.
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
    Eigen::MatrixXd M_;
    Eigen::VectorXd p_;
};

/**
 * \brief Control cost function class.
 * This cost function looks for a minimization of the control.
 * Mathematically, it is \f$(NU-p)^TW_U(NU-p) \Leftrightarrow sum_k w_u\|Nu_k-p\|^2\f$.
 */
class COPRA_DLLAPI ControlCost final : public CostFunction {
public:
    /**
     * \brief Constructor of the control cost function.
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
    Eigen::MatrixXd N_;
    Eigen::VectorXd p_;
};

/**
 * \brief Mixed cost function class.
 * This cost function looks for a minimization of a linear combination of trajectory and control.
 * Mathematically, it is \f$(MX+NU-p)^TW_M(MX+NU-p) \Leftrightarrow sum_k w_m\|Mx_k+Nu_k-p\|^2\f$.
 */
class COPRA_DLLAPI MixedCost final : public CostFunction {
public:
    /**
     * \brief Constructor of the mixed cost function.
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
    Eigen::MatrixXd M_, N_;
    Eigen::VectorXd p_;
};

} // namespace copra
