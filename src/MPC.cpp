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
#include "MPC.h"

// stl
#include <algorithm>
#include <exception>
#include <numeric>

// mpc
#include "PreviewSystem.h"
#include "constraints.h"
#include "costFunctions.h"

namespace mpc {

MPC::Constraints::Constraints()
    : nrEqConstr(0)
    , nrIneqConstr(0)
    , spConstr()
    , spEqConstr()
    , spIneqConstr()
    , spBoundConstr()
{
}

void MPC::Constraints::updateNr()
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

void MPC::Constraints::clear()
{
    nrEqConstr = 0;
    nrIneqConstr = 0;
    spConstr.clear();
    spEqConstr.clear();
    spIneqConstr.clear();
    spBoundConstr.clear();
}

/*************************************************************************************************
 *                                             MPC                                               *
 *************************************************************************************************/

MPC::MPC(SolverFlag sFlag)
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
    , solveTime_()
    , solveAndBuildTime_()
{
}

MPC::MPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag)
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
    , solveTime_()
    , solveAndBuildTime_()
{
    lb_.setConstant(-std::numeric_limits<double>::max());
    ub_.setConstant(std::numeric_limits<double>::max());
}

void MPC::selectQPSolver(SolverFlag flag)
{
    sol_ = solverFactory(flag);
}

void MPC::initializeController(const std::shared_ptr<PreviewSystem>& ps)
{
    ps_ = ps;
    clearConstraintMatrices();

    Q_.resize(ps_->fullUDim, ps_->fullUDim);
    c_.resize(ps_->fullUDim);
}

bool MPC::solve()
{
    using namespace std::chrono;
    auto sabTime = high_resolution_clock::now();

    updateSystem();
    makeQPForm();
    sol_->SI_problem(ps_->fullUDim, constraints_.nrEqConstr, constraints_.nrIneqConstr);

    auto sTime = high_resolution_clock::now();
    bool success = sol_->SI_solve(Q_, c_, Aeq_, beq_, Aineq_, bineq_, lb_, ub_);
    solveTime_ = duration_cast<duration<double> >(high_resolution_clock::now() - sTime);

    checkDeleteCostsAndConstraints();
    if (!success)
        sol_->SI_inform();

    solveAndBuildTime_ = duration_cast<duration<double> >(high_resolution_clock::now() - sabTime);
    return success;
}

const Eigen::VectorXd& MPC::control() const noexcept
{
    return sol_->SI_result();
}

Eigen::VectorXd MPC::trajectory() const noexcept
{
    return ps_->Phi * ps_->x0 + ps_->Psi * control() + ps_->xi;
}

double MPC::solveTime() const noexcept
{
    return solveTime_.count();
}

double MPC::solveAndBuildTime() const noexcept
{
    return solveAndBuildTime_.count();
}

void MPC::addCost(const std::shared_ptr<CostFunction>& costFun)
{
    costFun->initializeCost(*ps_);
    spCost_.emplace_back(costFun);
}

void MPC::addConstraint(const std::shared_ptr<Constraint>& constr)
{
    constr->initializeConstraint(*ps_);
    addConstraintByType(constr);
}

void MPC::resetConstraints() noexcept
{
    constraints_.clear();
    clearConstraintMatrices();
}

/*
 *  Protected methods
 */

void MPC::addConstraintByType(const std::shared_ptr<Constraint>& constr)
{
    switch (constr->constraintType()) {
    case ConstraintFlag::EqualityConstraint: {
        constraints_.nrEqConstr += constr->nrConstr();
        // DownCasting to std::shared_ptr<EqIneqConstraint>
        // This is a safe operation since we know that the object is a derived class of a EqIneqConstraint
        constraints_.spEqConstr.emplace_back(std::static_pointer_cast<EqIneqConstraint>(constr));
    ub_.setConstant(std::numeric_limits<double>::max());
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

void MPC::clearConstraintMatrices()
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

void MPC::updateSystem()
{
    // Reset the QP variables
    Q_.setIdentity();
    Q_ *= 1e-6; // Ensure that Q is positive definite if no cost functions are used
    c_.setZero();

    // Update the system
    if (!ps_->isUpdated)
        ps_->updateSystem();

    // Update the constraints
    constraints_.updateNr();
    Aeq_.resize(constraints_.nrEqConstr, ps_->fullUDim);
    beq_.resize(constraints_.nrEqConstr);
    Aineq_.resize(constraints_.nrIneqConstr, ps_->fullUDim);
    bineq_.resize(constraints_.nrIneqConstr);

    for (auto& cstr : constraints_.spConstr)
        cstr->update(*ps_);

    // Update the costs
    for (auto& sp : spCost_)
        sp->update(*ps_);
}

void MPC::makeQPForm()
{
    for (auto& cost : spCost_) {
        Q_ += cost->Q();
        c_ += cost->c();
    }

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

void MPC::checkDeleteCostsAndConstraints()
{
    auto check = [](auto& sp, bool useWarn = false, int delLimit = 2) {
        for (auto it = sp.begin(); it != sp.end();) {
            if ((*it).use_count() <= delLimit) {
                CONSTRAINT_DELETION_WARN(useWarn, "%s%s%s", "A '", (*it)->name().c_str(),
                    "' has been destroyed.\nIt has been removed from the controller\n");
                (void)useWarn;
                it = sp.erase(it);
            } else {
                ++it;
            }
        }
    };

    check(constraints_.spConstr, true);
    check(constraints_.spEqConstr);
    check(constraints_.spIneqConstr);
    check(constraints_.spBoundConstr);
    check(constraints_.spBoundConstr);
    check(spCost_, true, 1);
}

} // namespace mpc