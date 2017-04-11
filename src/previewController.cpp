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
#include "previewController.h"

// stl
#include <algorithm>
#include <exception>
#include <numeric>

// mpc
#include "PreviewSystem.h"
#include "constraints.h"
#include "debugUtils.h"

namespace mpc {

MPCTypeFull::Constraints::Constraints()
    : nrEqConstr(0)
    , nrIneqConstr(0)
    , spConstr()
    , spEqConstr()
    , spIneqConstr()
    , spBoundConstr()
{
}

void MPCTypeFull::Constraints::updateNr()
{
    nrEqConstr = 0;
    nrIneqConstr = 0;

    auto countConstr = [](auto& spc, int& nreq) {
        for (auto& sp : spc)
            nreq += sp->nrConstr();
    };

    countConstr(spEqConstr, nrEqConstr);
    countConstr(spIneqConstr, nrIneqConstr);
}

void MPCTypeFull::Constraints::clear()
{
    nrEqConstr = 0;
    nrIneqConstr = 0;
    spConstr.clear();
    spEqConstr.clear();
    spIneqConstr.clear();
    spBoundConstr.clear();
}

/*************************************************************************************************
 *                                         MPCTypeFull                                           *
 *************************************************************************************************/

MPCTypeFull::MPCTypeFull(SolverFlag sFlag)
    : ps_(nullptr)
    , sol_(solverFactory(sFlag))
    , constraints_()
    , Q_()
    , Aineq_()
    , Aeq_()
    , c_()
    , bineq_()
    , beq_()
    , lb_()
    , ub_()
    , Wx_()
    , Wu_()
    , solveTime_()
    , solveAndBuildTime_()
{
}

MPCTypeFull::MPCTypeFull(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag)
    : ps_(ps)
    , sol_(solverFactory(sFlag))
    , constraints_()
    , Q_(ps_->fullUDim, ps_->fullUDim)
    , Aineq_(0, ps_->fullUDim)
    , Aeq_(0, ps_->fullUDim)
    , c_(ps_->fullUDim)
    , bineq_()
    , beq_()
    , lb_(ps_->fullUDim)
    , ub_(ps_->fullUDim)
    , Wx_(ps_->fullXDim)
    , Wu_(ps_->fullUDim)
    , solveTime_()
    , solveAndBuildTime_()
{
    Wx_.setOnes();
    Wu_.setOnes();
    lb_.setConstant(-std::numeric_limits<double>::max());
    ub_.setConstant(std::numeric_limits<double>::max());
    bineq_.resize(0);
    beq_.resize(0);
}

void MPCTypeFull::selectQPSolver(SolverFlag flag)
{
    sol_ = solverFactory(flag);
}

void MPCTypeFull::initializeController(const std::shared_ptr<PreviewSystem>& ps)
{
    ps_ = ps;
    clearConstraintMatrices();

    Q_.resize(ps_->fullUDim, ps_->fullUDim);
    c_.resize(ps_->fullUDim);
    Wx_.resize(ps_->fullXDim);
    Wu_.resize(ps_->fullUDim);
    Wx_.setOnes();
    Wu_.setOnes();
}

bool MPCTypeFull::solve()
{
    solveAndBuildTime_.start();
    updateSystem();
    makeQPForm();
    sol_->SI_problem(ps_->fullUDim, constraints_.nrEqConstr, constraints_.nrIneqConstr);
    solveTime_.start();
    bool success = sol_->SI_solve(Q_, c_, Aeq_, beq_, Aineq_, bineq_, lb_, ub_);
    solveTime_.stop();
    checkDeleteConstraints();
    solveAndBuildTime_.stop();
    if (!success)
        sol_->SI_inform();

    return success;
}

const Eigen::VectorXd& MPCTypeFull::control() const noexcept
{
    return sol_->SI_result();
}

Eigen::VectorXd MPCTypeFull::trajectory() const noexcept
{
    return ps_->Phi * ps_->x0 + ps_->Psi * control() + ps_->xi;
}

boost::timer::cpu_times MPCTypeFull::solveTime() const noexcept
{
    return solveTime_.elapsed();
}

boost::timer::cpu_times MPCTypeFull::solveAndBuildTime() const noexcept
{
    return solveAndBuildTime_.elapsed();
}

void MPCTypeFull::weights(const Eigen::VectorXd& Wx, const Eigen::VectorXd& Wu)
{
    if (Wx.rows() == ps_->xDim)
        for (auto i = 0; i < ps_->nrXStep; ++i)
            Wx_.segment(i * ps_->xDim, ps_->xDim) = Wx;
    else if (Wx.rows() == ps_->fullXDim)
        Wx_ = Wx;
    else
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSXDim("Wx", Wx, ps_.get()));

    if (Wu.rows() == ps_->uDim)
        for (auto i = 0; i < ps_->nrUStep; ++i)
            Wu_.segment(i * ps_->uDim, ps_->uDim) = Wu;
    else if (Wu.rows() == ps_->fullUDim)
        Wu_ = Wu;
    else
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSUDim("Wu", Wu, ps_.get()));
}

void MPCTypeFull::weights(double Wx, double Wu)
{
    Wx_.setConstant(Wx);
    Wu_.setConstant(Wu);
}

void MPCTypeFull::addConstraint(const std::shared_ptr<Constraint>& constr)
{
    constr->initializeConstraint(*ps_);
    addConstraintByType(constr);
}

void MPCTypeFull::resetConstraints() noexcept
{
    constraints_.clear();
    clearConstraintMatrices();
}

/*
 *  Protected methods
 */

void MPCTypeFull::addConstraintByType(const std::shared_ptr<Constraint>& constr)
{
    switch (constr->constraintType()) {
    case ConstraintFlag::EqualityConstraint: {
        constraints_.nrEqConstr += constr->nrConstr();
        // DownCasting to std::shared_ptr<EqIneqConstraint>
        // This is a safe operation since we know that the object is a derived class of a EqIneqConstraint
        constraints_.spEqConstr.emplace_back(std::static_pointer_cast<EqIneqConstraint>(constr));
    } break;
    case ConstraintFlag::InequalityConstraint: {
        constraints_.nrIneqConstr += constr->nrConstr();
        // DownCasting to std::shared_ptr<EqIneqConstraint>
        // This is a safe operation since we know that the object is a derived class of a EqIneqConstraint
        constraints_.spIneqConstr.emplace_back(std::static_pointer_cast<EqIneqConstraint>(constr));
    } break;
    case ConstraintFlag::BoundConstraint: {
        // DownCasting to std::shared_ptr<ControlBoundConstraint>
        // This is a safe operation since we know that the object is a ControlBoundConstraint
        constraints_.spBoundConstr.emplace_back(std::static_pointer_cast<ControlBoundConstraint>(constr));
    } break;
    default:
        return;
    }
    constraints_.spConstr.emplace_back(constr);
}

void MPCTypeFull::clearConstraintMatrices()
{
    assert(ps_ != nullptr);

    Aineq_.resize(0, ps_->fullUDim);
    Aeq_.resize(0, ps_->fullUDim);
    bineq_.resize(0);
    beq_.resize(0);
    lb_.resize(ps_->fullUDim);
    ub_.resize(ps_->fullUDim);
    lb_.setConstant(-std::numeric_limits<double>::max());
    ub_.setConstant(std::numeric_limits<double>::max());
}

void MPCTypeFull::updateSystem()
{
    ps_->updateSystem();

    constraints_.updateNr();
    Aeq_.resize(constraints_.nrEqConstr, ps_->fullUDim);
    beq_.resize(constraints_.nrEqConstr);
    Aineq_.resize(constraints_.nrIneqConstr, ps_->fullUDim);
    bineq_.resize(constraints_.nrIneqConstr);
    lb_.resize(ps_->fullUDim);
    ub_.resize(ps_->fullUDim);
    lb_.setConstant(-std::numeric_limits<double>::max());
    ub_.setConstant(std::numeric_limits<double>::max());

    for (auto& sp : constraints_.spConstr)
        sp->update(*ps_);
}

void MPCTypeFull::makeQPForm()
{
    //REDO
    Q_.noalias() = ps_->Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * ps_->Psi + Eigen::MatrixXd(Wu_.asDiagonal());
    c_.noalias() = ps_->Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (ps_->Phi * ps_->x0 - ps_->xd + ps_->xi);

    int nrLines = 0;
    // Get Equality constraints
    for (auto& cstr : constraints_.spEqConstr) {
        Aeq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        beq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }

    // Get Inequality constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spIneqConstr) {
        Aineq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        bineq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }

    // Get Bound constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spBoundConstr) {
        lb_.segment(nrLines, cstr->nrConstr()) = cstr->lower();
        ub_.segment(nrLines, cstr->nrConstr()) = cstr->upper();
        nrLines += cstr->nrConstr();
    }
}

void MPCTypeFull::checkDeleteConstraints()
{
    auto checkConstr = [](auto& spc, bool useWarn = false) {
        for (auto it = spc.begin(); it != spc.end();) {
            if ((*it).use_count() <= 2) {
                CONSTRAINT_DELETION_WARN(useWarn, "%s%s%s", "A '", (*it)->name().c_str(),
                    "' has been destroyed.\nThe constraint has been removed from the controller\n");
                (void)useWarn;
                it = spc.erase(it);
            } else {
                ++it;
            }
        }
    };

    checkConstr(constraints_.spConstr, true);
    checkConstr(constraints_.spEqConstr);
    checkConstr(constraints_.spIneqConstr);
    checkConstr(constraints_.spBoundConstr);
}

/*************************************************************************************************
 *                                         MPCTypeLast                                           *
 *************************************************************************************************/

MPCTypeLast::MPCTypeLast(SolverFlag sFlag)
    : MPCTypeFull::MPCTypeFull(sFlag)
{
}

MPCTypeLast::MPCTypeLast(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag)
    : MPCTypeFull::MPCTypeFull(ps, sFlag)
{
    Wx_.resize(ps_->xDim);
    Wx_.setOnes();
}

void MPCTypeLast::initializeController(const std::shared_ptr<PreviewSystem>& ps)
{
    ps_ = ps;
    Q_.resize(ps_->fullUDim, ps_->fullUDim);
    Aineq_.resize(0, ps_->fullUDim);
    Aeq_.resize(0, ps_->fullUDim);
    c_.resize(ps_->fullUDim);
    bineq_.resize(0);
    beq_.resize(0);
    lb_.resize(ps_->fullUDim);
    ub_.resize(ps_->fullUDim);
    lb_.setConstant(-std::numeric_limits<double>::max());
    ub_.setConstant(std::numeric_limits<double>::max());
    Wx_.resize(ps_->xDim);
    Wu_.resize(ps_->fullUDim);
    Wx_.setOnes();
    Wu_.setOnes();
}

void MPCTypeLast::weights(const Eigen::VectorXd& Wx, const Eigen::VectorXd& Wu)
{
    if (Wx.rows() == ps_->xDim)
        Wx_ = Wx;
    else
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnDim("Wu", Wu, ps_->xDim));
        
    if (Wu.rows() == ps_->uDim)
        for (auto i = 0; i < ps_->nrUStep; ++i)
            Wu_.segment(i * ps_->uDim, ps_->uDim) = Wu;
    else if (Wu.rows() == ps_->fullUDim)
        Wu_ = Wu;
    else
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSUDim("Wu", Wu, ps_.get()));
}

/*
 *  Protected methods
 */

void MPCTypeLast::makeQPForm()
{
    auto xDim = ps_->xDim;
    const Eigen::MatrixXd& psi = ps_->Psi.bottomRows(xDim);
    Q_.noalias() = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * psi + Eigen::MatrixXd(Wu_.asDiagonal());
    c_.noalias() = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (ps_->Phi.bottomRows(xDim) * ps_->x0 - ps_->xd.tail(xDim) + ps_->xi.tail(xDim));

    int nrLines = 0;
    // Get Equality constraints
    for (auto& cstr : constraints_.spEqConstr) {
        Aeq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        beq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }

    // Get Inequality constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spIneqConstr) {
        Aineq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        bineq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }

    // Get Bound constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spBoundConstr) {
        lb_.segment(nrLines, cstr->nrConstr()) = cstr->lower();
        ub_.segment(nrLines, cstr->nrConstr()) = cstr->upper();
        nrLines += cstr->nrConstr();
    }
}

} // namespace pc