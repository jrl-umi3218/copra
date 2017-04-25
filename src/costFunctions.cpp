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

namespace mpc {

/*************************************************************************************************
 *                                        Cost Function                                          *
 *************************************************************************************************/

CostFunction::CostFunction(std::string&& name)
    : name_(std::move(name))
{
}

void CostFunction::autoSpan()
{
}

void CostFunction::initializeCost(const std::shared_ptr<PreviewSystem> ps)
{
    Q_.resize(ps->fullUdim, ps->fullUDim);
    c_.resize(ps->fullUdim);
}

/*************************************************************************************************
 *                                  Trajectory Cost Function                                     *
 *************************************************************************************************/

void TrajectoryCost::initializeCost(const PreviewSystem& ps)
{
    using CostFunction::initializeConstraint(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));

    if (M_.cols() == ps.xDim) {
        AutoSpan::spanMatrix(M_, M_.rows() * ps.nrXStep); // Use autospan but this will change
        AutoSpan::spanMatrix(weights_, weights_.rows() * ps.nrXStep); // Use autospan but this will change
        AutoSpan::spanVector(p_, p_.rows() * ps.nrXStep); // Use autospan but this will change
    } else if (M_.cols() == ps.fullXDim) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXDim("M", M_, &ps));
    }
}

void TrajectoryCost::update(const PreviewSystem& ps)
{
    Eigen::MatrixXd tmp = M_ * ps.Psi;
    Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
    c_.noalias() = (M_ * (ps.Phi * ps.x0 + ps.xi) + p_).transpose() * weights_ * tmp;
}

/*************************************************************************************************
 *                                    Target Cost Function                                       *
 *************************************************************************************************/

void TargetCost::initializeCost(const PreviewSystem& ps)
{
    using CostFunction::initializeConstraint(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));

    if (M_.cols() != ps.xDim)
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSxDim("M", M_, &ps));
}

void TargetCost::update(const PreviewSystem& ps)
{
    Eigen::MatrixXd tmp = M_ * ps.Psi.bottomRows(ps.xDim);
    Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
    c_.noalias() = (M_ * (ps.Phi.bottomRows(ps.xDim) * ps.x0 + ps.xi.bottomRows(ps.xDim)) + p_).transpose() * weights_ * tmp;
}

/*************************************************************************************************
 *                                   Control Cost Function                                       *
 *************************************************************************************************/

void ControlCost::initializeCost(const PreviewSystem& ps)
{
    using CostFunction::initializeConstraint(ps);
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
    if (fullSize_) {
        Q_.noalias() = N_.transpose() * weights_.asDiagonal() * N_;
        c_.noalias() = p_.transpose() * weights_.asDiagonal() * N_;
    } else {
        Eigen::MatrixXd mat = N_.transpose() * weights_.asDiagonal() * N_;
        Eigen::MatrixXd vec = p_.transpose() * weights_.asDiagonal() * N_;
        for (int i = 0; i < ps.nbUStep; ++i) {
            Q_.block(i * ps.uDim, i * ps.uDim, ps.uDim, ps.uDim) = mat;
            c_.segment(i * ps.uDim, ps.uDim) = vec;
        }
    }
}

/*************************************************************************************************
 *                               Mixed Trajectory Cost Function                                  *
 *************************************************************************************************/

void MixedTrajectoryCost::initializeCost(const PreviewSystem& ps)
{
    using CostFunction::initializeConstraint(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));
    if (N_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("N", "p", N_, p_));

    if (M_.cols() == ps.xDim && N_.cols() == ps.uDim) {
        AutoSpan::spanMatrix(M_, M_.rows() * ps.nrUStep, 1); // Use autospan but this will change
        AutoSpan::spanMatrix(N_, N_.rows() * ps.nrUStep); // Use autospan but this will change
        AutoSpan::spanMatrix(weights_, weights_.rows() * ps.nrUStep); // Use autospan but this will change
        AutoSpan::spanVector(p_, p_.rows() * ps.nrUStep); // Use autospan but this will change
    } else if (M_.cols() != ps.fullXDim || N_.cols() != ps.fullUDim) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXUDim("M", "N", M_, N_, &ps));
    }
}

void MixedTrajectoryCost::update(const PreviewSystem& ps)
{
    Eigen::MatrixXd tmp = M_ * ps.Psi + N_;
    Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
    c_.noalias() = (M_ * (ps.Phi * ps.x0 + ps.xi) + p_).transpose() * weights_.asDiagonal() * tmp;
}

/*************************************************************************************************
 *                                 Mixed Target Cost Function                                    *
 *************************************************************************************************/

void MixedTargetCost::initializeCost(const PreviewSystem& ps)
{
    using CostFunction::initializeConstraint(ps);
    if (M_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("M", "p", M_, p_));
    if (N_.rows() != p_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsAskAutoSpan("N", "p", N_, p_));

    if (M_.cols() != ps.xDim || N_.cols() != ps.uDim)
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSxuDim("M", "N", M_, N_, &ps));
}

void MixedTargetCost::update(const PreviewSystem& ps)
{
    Eigen::MatrixXd tmp = M_ * ps.Psi.bottomRows(ps.xDim) + N_;
    Q_.noalias() = tmp.transpose() * weights_.asDiagonal() * tmp;
    c_.noalias() = (M_ * (ps.Phi.bottomRows(ps.xDim) * ps.x0 + ps.xi.bottomRows(ps.xDim)) + p_).transpose() * weights_.asDiagonal() * tmp;
}

} // namespace mpc