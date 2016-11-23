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

#include "PreviewController.h"

#include "Constraints.h"
#include "PreviewSystem.h"
#include "debugUtils.h"
#include <algorithm>
#include <exception>
#include <numeric>
#include <sstream>

namespace mpc {

/*************************************************************************************************
 *                                         MPCTypeFull *
 *************************************************************************************************/

MPCTypeFull::MPCTypeFull(SolverFlag sFlag)
    : nrConstr_(0)
    , ps_()
    , wpc_()
    , sol_(solverFactory(sFlag))
    , Q_()
    , AInEq_()
    , c_()
    , bInEq_()
    , Wx_()
    , Wu_()
    , solveTime_()
    , solveAndBuildTime_()
{
}

MPCTypeFull::MPCTypeFull(const std::shared_ptr<PreviewSystem>& ps,
    SolverFlag sFlag)
    : nrConstr_(0)
    , ps_(ps)
    , wpc_()
    , sol_(solverFactory(sFlag))
    , Q_(ps_->fullUDim, ps_->fullUDim)
    , AInEq_(3 * ps_->fullUDim, 3 * ps_->fullUDim)
    , // max size is 1 equality (= 2 inequalities) + 1 inequality
    c_(ps_->fullUDim)
    , bInEq_(3 * ps_->fullUDim)
    , Wx_(ps_->fullXDim)
    , Wu_(ps_->fullUDim)
    , solveTime_()
    , solveAndBuildTime_()
{
    Wx_.setOnes();
    Wu_.setOnes();
}

void MPCTypeFull::selectQPSolver(SolverFlag flag)
{
    sol_ = solverFactory(flag);
}

void MPCTypeFull::initializeController(const std::shared_ptr<PreviewSystem>& ps)
{
    ps_ = ps;
    Q_.resize(ps_->fullUDim, ps_->fullUDim);
    AInEq_.resize(3 * ps_->fullUDim, 3 * ps_->fullUDim); // max size is 1 equality
    // (= 2 inequalities) + 1
    // inequality
    c_.resize(ps_->fullUDim);
    bInEq_.resize(3 * ps_->fullUDim);
    Wx_.resize(ps_->fullXDim);
    Wu_.resize(ps_->fullUDim);
    Wx_.setOnes();
    Wu_.setOnes();
}

bool MPCTypeFull::solve()
{
    solveAndBuildTime_.start();
    checkConstraints();
    updateSystem();
    makeQPForm();
    sol_->SI_problem(ps_->fullUDim, 0, nrConstr_);
    solveTime_.start();
    bool success = sol_->SI_solve(
        Q_, c_, Eigen::MatrixXd::Zero(0, ps_->fullUDim), Eigen::VectorXd::Zero(0),
        AInEq_, bInEq_, Eigen::VectorXd::Constant(ps_->fullUDim, -std::numeric_limits<double>::infinity()),
        Eigen::VectorXd::Constant(ps_->fullUDim, std::numeric_limits<double>::infinity()));
    solveTime_.stop();
    solveAndBuildTime_.stop();
    if (!success)
        sol_->SI_inform();

    return success;
}

const Eigen::VectorXd&
MPCTypeFull::control() const noexcept
{
    return sol_->SI_result();
}

Eigen::VectorXd
MPCTypeFull::trajectory() const noexcept
{
    return ps_->Phi * ps_->x0 + ps_->Psi * control() + ps_->xi;
}

boost::timer::cpu_times
MPCTypeFull::solveTime() const noexcept
{
    return solveTime_.elapsed();
}

boost::timer::cpu_times
MPCTypeFull::solveAndBuildTime() const noexcept
{
    return solveAndBuildTime_.elapsed();
}

void MPCTypeFull::weights(const Eigen::VectorXd& Wx, const Eigen::VectorXd& Wu)
{
    assert(Wx.rows() == ps_->xDim);
    assert(Wu.rows() == ps_->uDim);

    for (auto i = 0; i < ps_->nrStep; ++i) {
        Wx_.segment(i * ps_->xDim, ps_->xDim) = Wx;
        Wu_.segment(i * ps_->uDim, ps_->uDim) = Wu;
    }
}

void MPCTypeFull::addConstraint(std::shared_ptr<Constraint> constr)
{
    wpc_.emplace_back(std::move(constr), constr->name());
    constr->initializeConstraint(*ps_);
    nrConstr_ += constr->nrConstr();

    AInEq_.resize(nrConstr_, ps_->fullUDim);
    bInEq_.resize(nrConstr_);
}

void MPCTypeFull::resetConstraints() noexcept
{
    nrConstr_ = 0;
    wpc_.clear();
}

/*
 *  Protected methods
 */

void MPCTypeFull::updateSystem()
{
    auto xDim = ps_->xDim;
    auto uDim = ps_->uDim;

    ps_->Phi.block(0, 0, xDim, xDim) = ps_->A;
    ps_->Psi.block(0, 0, xDim, uDim) = ps_->B;
    ps_->xi.segment(0, xDim) = ps_->d;

    for (auto i = 1; i < ps_->nrStep; ++i) {
        ps_->Phi.block(i * xDim, 0, xDim, xDim) = ps_->A * ps_->Phi.block((i - 1) * xDim, 0, xDim, xDim);
        for (auto j = 0; j < i; ++j) {
            ps_->Psi.block(i * xDim, j * uDim, xDim, uDim) = ps_->A * ps_->Psi.block((i - 1) * xDim, j * uDim, xDim, uDim);
        }
        ps_->Psi.block(i * xDim, i * uDim, xDim, uDim) = ps_->B;
        ps_->xi.segment(i * xDim, xDim) = ps_->A * ps_->xi.segment((i - 1) * xDim, xDim) + ps_->d;
    }

    for (auto& wpc : wpc_)
        wpc.first.lock()->update(*ps_);
}

void MPCTypeFull::makeQPForm()
{
    Q_ = ps_->Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * ps_->Psi + Eigen::MatrixXd(Wu_.asDiagonal());
    c_ = ps_->Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (ps_->Phi * ps_->x0 - ps_->xd + ps_->xi);

    int nrLines = 0;
    for (auto& wpc : wpc_) {
        auto cstr = wpc.first.lock();
        AInEq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        bInEq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }
}

void MPCTypeFull::checkConstraints()
{
    bool needResizing = false;
    for (auto itr = wpc_.begin(); itr != wpc_.end();) {
        if ((*itr).first.expired()) {
            std::stringstream ss;
            ss << "Dangling pointer to constrain.\nA '" << (*itr).second
               << "' has been destroyed unexpectedly.\n"
                  "The constraint has been removed from the controller";
            DEBUG_WARN(ss.str());
            itr = wpc_.erase(itr);
            needResizing = true;
        } else {
            ++itr;
        }
    }

    if (needResizing) {
        nrConstr_ = 0;
        for (auto& wpc : wpc_)
            nrConstr_ += wpc.first.lock()->nrConstr();
        AInEq_.resize(nrConstr_, ps_->fullUDim);
        bInEq_.resize(nrConstr_);
    }
}

/*************************************************************************************************
 *                                         MPCTypeLast *
 *************************************************************************************************/

MPCTypeLast::MPCTypeLast(SolverFlag sFlag)
    : MPCTypeFull::MPCTypeFull(sFlag)
{
}

MPCTypeLast::MPCTypeLast(const std::shared_ptr<PreviewSystem>& ps,
    SolverFlag sFlag)
    : MPCTypeFull::MPCTypeFull(ps, sFlag)
{
    Wx_.resize(ps_->xDim);
    Wx_.setOnes();
}

void MPCTypeLast::initializeController(const std::shared_ptr<PreviewSystem>& ps)
{
    ps_ = ps;
    Q_.resize(ps_->fullUDim, ps_->fullUDim);
    AInEq_.resize(3 * ps_->fullUDim, 3 * ps_->fullUDim); // max size is 1 equality
    // (= 2 inequalities) + 1
    // inequality
    c_.resize(ps_->fullUDim);
    bInEq_.resize(3 * ps_->fullUDim);
    Wx_.resize(ps_->xDim);
    Wu_.resize(ps_->fullUDim);
    Wx_.setOnes();
    Wu_.setOnes();
}

void MPCTypeLast::weights(const Eigen::VectorXd& Wx, const Eigen::VectorXd& Wu)
{
    assert(Wx.rows() == ps_->xDim);
    assert(Wu.rows() == ps_->uDim);

    Wx_ = Wx;
    for (auto i = 0; i < ps_->nrStep; ++i)
        Wu_.segment(i * ps_->uDim, ps_->uDim) = Wu;
}

/*
 *  Protected methods
 */

void MPCTypeLast::makeQPForm()
{
    auto xDim = ps_->xDim;
    const Eigen::MatrixXd& psi = ps_->Psi.bottomRows(xDim);
    Q_ = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * psi + Eigen::MatrixXd(Wu_.asDiagonal());
    c_ = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (ps_->Phi.bottomRows(xDim) * ps_->x0 - ps_->xd.tail(xDim) + ps_->xi.tail(xDim));

    int nrLines = 0;
    checkConstraints();
    for (auto& wpc : wpc_) {
        auto cstr = wpc.first.lock();
        AInEq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        bInEq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }
}

} // namespace pc