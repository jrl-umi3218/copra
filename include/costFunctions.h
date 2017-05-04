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

// Eigen
#include <Eigen/Core>

// mpc
#include "config.hh"
#include "debugUtils.h"
#include "typedefs.h"

namespace mpc {

class MPC_DLLAPI CostFunction {
public:
    CostFunction(std::string&& name);
    virtual ~CostFunction() = default;

    virtual void autoSpan();
    virtual void update(const PreviewSystem& ps) = 0;
    virtual void initializeCost(const PreviewSystem& ps);

    /**
     * Set the weights of the system.
     * Perform a move semantic if a fullsize vector is given as rvalue (this is faster).
     * \param Weights of to set.
     * \throw Throw a std::domain_error if weights is badly dimension.
     */
    template <typename TVec, typename = std::enable_if_t<!is_all_arithmetic<TVec>::value> >
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
    template <typename T, typename = std::enable_if_t<is_all_arithmetic<T>::value> >
    void weights(T weight)
    {
        weights_.setConstant(weight);
    }

    const std::string& name() const noexcept
    {
        return name_;
    }

    const Eigen::MatrixXd& Q() const noexcept
    {
        return Q_;
    }

    const Eigen::VectorXd& c() const noexcept
    {
        return c_;
    }

protected:
    std::string name_;
    bool fullSizeEntry_;
    Eigen::MatrixXd Q_;
    Eigen::VectorXd c_;
    Eigen::VectorXd weights_;
};

class TrajectoryCost final : public CostFunction {
public:
    template <typename TMat, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value> >
    TrajectoryCost(TMat&& M, TVec&& p)
        : CostFunction("TrajectoryCost")
        , M_(std::forward<TMat>(M))
        , p_(std::forward<TVec>(p))
    {
        weights_ = Eigen::VectorXd::Ones(p_.rows());
    }

    void autoSpan() override;
    void update(const PreviewSystem& ps) override;
    void initializeCost(const PreviewSystem& ps) override;

private:
    Eigen::MatrixXd M_;
    Eigen::VectorXd p_;
};

class TargetCost final : public CostFunction {
public:
    template <typename TMat, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value> >
    TargetCost(TMat&& M, TVec&& p)
        : CostFunction("TargetCost")
        , M_(std::forward<TMat>(M))
        , p_(std::forward<TVec>(p))
    {
        weights_ = Eigen::VectorXd::Ones(p_.rows());
    }

    void update(const PreviewSystem& ps) override;
    void initializeCost(const PreviewSystem& ps) override;

private:
    Eigen::MatrixXd M_;
    Eigen::VectorXd p_;
};

class ControlCost final : public CostFunction {
public:
    template <typename TMat, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat, TVec>::value> >
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

class MixedCost final : public CostFunction {
public:
    template <typename TMat1, typename TMat2, typename TVec,
        typename = std::enable_if_t<!is_all_arithmetic<TMat1, TMat2, TVec>::value> >
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

} // namespace mpc