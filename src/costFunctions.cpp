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

// header
#include "costFunctions.h"

// mpc
#include "AutoSpan.h"

namespace mpc {

/*************************************************************************************************
 *                                        Cost Function                                          *
 *************************************************************************************************/

CostFunction::CostFunction(std::string&& name)
    : name_(std::move(name))
    , fullSizeEntry_(false)
{
}

void CostFunction::autoSpan()
{
}

void CostFunction::initializeCost(const PreviewSystem& ps)
{
    Q_.resize(ps.fullUDim, ps.fullUDim);
    c_.resize(ps.fullUDim);
}

/*************************************************************************************************
 *                                  Trajectory Cost Function                                     *
 *************************************************************************************************/

void TrajectoryCost::autoSpan()
{
    auto max_dim = std::max(M_.rows(), std::max(weights_.rows(), p_.rows()));
    AutoSpan::spanMatrix(M_, max_dim);
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void TrajectoryCost::initializeCost(const PreviewSystem& ps)
{
    CostFunction::initializeCost(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));

    if (M_.cols() == ps.xDim) {
        Q_.setZero();
        c_.setZero();
    } else if (M_.cols() == ps.fullXDim) {
        fullSizeEntry_ = true;
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXDim("M", M_, &ps));
    }
}

void TrajectoryCost::update(const PreviewSystem& ps)
{
    if (fullSizeEntry_) {
        Eigen::MatrixXd tmp {M_ * ps.Psi};
        Q_ = tmp.transpose() * weights_.asDiagonal() * tmp;
        c_ = (M_ * (ps.Phi * ps.x0 + ps.xi) + p_).transpose() * weights_.asDiagonal() * tmp;
    } else {
        for (int i = 0; i < ps.nrXStep; ++i) { // Can be optimized. Lot of sums of zero here
            Eigen::MatrixXd tmp {M_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim)};
            Q_.noalias() += tmp.transpose() * weights_.asDiagonal() * tmp;
            c_.noalias() += (M_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim)) + p_).transpose() * weights_.asDiagonal() * tmp;
        }
    }
}

/*************************************************************************************************
 *                                    Target Cost Function                                       *
 *************************************************************************************************/

void TargetCost::initializeCost(const PreviewSystem& ps)
{
    CostFunction::initializeCost(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));

    if (M_.cols() != ps.xDim)
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSxDim("M", M_, &ps));
}

void TargetCost::update(const PreviewSystem& ps)
{
    Eigen::MatrixXd tmp {M_ * ps.Psi.bottomRows(ps.xDim)};
    Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
    c_.noalias() = (M_ * (ps.Phi.bottomRows(ps.xDim) * ps.x0 + ps.xi.bottomRows(ps.xDim)) + p_).transpose() * weights_.asDiagonal() * tmp;
}

/*************************************************************************************************
 *                                   Control Cost Function                                       *
 *************************************************************************************************/

void ControlCost::autoSpan()
{
    auto max_dim = std::max(N_.rows(), std::max(weights_.rows(), p_.rows()));
    AutoSpan::spanMatrix(N_, max_dim);
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void ControlCost::initializeCost(const PreviewSystem& ps)
{
    CostFunction::initializeCost(ps);
    if (N_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("N", "p", N_, p_));

    if (N_.cols() == ps.uDim)
        Q_.setZero();
    else if (N_.cols() == ps.fullUDim)
        fullSizeEntry_ = true;
    else
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSUDim("N", N_, &ps));
}

void ControlCost::update(const PreviewSystem& ps)
{
    if (fullSizeEntry_) {
        Q_.noalias() = N_.transpose() * weights_.asDiagonal() * N_;
        c_.noalias() = p_.transpose() * weights_.asDiagonal() * N_;
    } else {
        Eigen::MatrixXd mat {N_.transpose() * weights_.asDiagonal() * N_};
        Eigen::VectorXd vec {p_.transpose() * weights_.asDiagonal() * N_};
        for (int i = 0; i < ps.nrUStep; ++i) {
            Q_.block(i * ps.uDim, i * ps.uDim, ps.uDim, ps.uDim) = mat;
            c_.segment(i * ps.uDim, ps.uDim) = vec;
        }
    }
}

/*************************************************************************************************
 *                                    Mixed Cost Function                                        *
 *************************************************************************************************/

void MixedCost::autoSpan()
{
    auto max_dim = std::max(M_.rows(), std::max(N_.rows(), std::max(weights_.rows(), p_.rows())));
    AutoSpan::spanMatrix(M_, max_dim, 1); // This is tricky. Has X and U are not the same dimensions, we need to handle it.
    AutoSpan::spanMatrix(N_, max_dim);
    AutoSpan::spanVector(p_, max_dim);
    AutoSpan::spanVector(weights_, max_dim);
}

void MixedCost::initializeCost(const PreviewSystem& ps)
{
    CostFunction::initializeCost(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));
    if (N_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("N", "p", N_, p_));

    if (M_.cols() == ps.xDim && N_.cols() == ps.uDim) {
        Q_.setZero();
        c_.setZero();
    } else if (M_.cols() == ps.fullXDim && N_.cols() == ps.fullUDim) {
        fullSizeEntry_ = true;
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXUDim("M", "N", M_, N_, &ps));
    }
}

void MixedCost::update(const PreviewSystem& ps)
{
    if (fullSizeEntry_) {
        Eigen::MatrixXd tmp = M_ * ps.Psi + N_;
        Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
        c_.noalias() = (M_ * (ps.Phi * ps.x0 + ps.xi) + p_).transpose() * weights_.asDiagonal() * tmp;
    } else {
        for (int i = 0; i < ps.nrUStep; ++i) { // Can be optimized. Lot of sums of zero here
            Eigen::MatrixXd tmp {M_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim)};
            tmp.block(0, i * ps.uDim, M_.rows(), ps.uDim).noalias() += N_;
            Q_.noalias() += tmp.transpose() * weights_.asDiagonal() * tmp;
            c_.noalias() += (M_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim)) + p_).transpose() * weights_.asDiagonal() * tmp;
        }
    }
}

} // namespace mpc
