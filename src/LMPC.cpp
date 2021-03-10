/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "LMPC.h"
#include "PreviewSystem.h"
#include "constraints.h"
#include "costFunctions.h"
#include <algorithm>
#include <exception>
#include <numeric>

namespace copra {

LMPC::Constraints::Constraints()
    : nrEqConstr(0)
    , nrIneqConstr(0)
{
}

void LMPC::Constraints::updateNr()
{
    nrEqConstr = 0;
    nrIneqConstr = 0;

    auto countConstr = [](auto& spc, int& nreq) {
        for (auto& sp : spc) {
            nreq += sp->nrConstr();
        }
    };

    countConstr(spEqConstr, nrEqConstr);
    countConstr(spIneqConstr, nrIneqConstr);
}

void LMPC::Constraints::clear()
{
    nrEqConstr = 0;
    nrIneqConstr = 0;
    spConstr.clear();
    spEqConstr.clear();
    spIneqConstr.clear();
    spBoundConstr.clear();
}

/*************************************************************************************************
 *                                             LMPC                                              *
 *************************************************************************************************/

LMPC::LMPC(SolverFlag sFlag)
    : ps_(nullptr)
    , sol_(solverFactory(sFlag))
{
}

LMPC::LMPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag)
    : sol_(solverFactory(sFlag))
{
    initializeController(ps);
}

void LMPC::selectQPSolver(SolverFlag flag)
{
    sol_ = solverFactory(flag);
}

void LMPC::useSolver(std::unique_ptr<SolverInterface>&& solver)
{
    sol_ = std::move(solver);
}

void LMPC::initializeController(const std::shared_ptr<PreviewSystem>& ps)
{
    ps_ = ps;
    clearConstraintMatrices();
    updateQPMatrixSize();
}

bool LMPC::solve()
{
    using namespace std::chrono;
    auto sabTime = high_resolution_clock::now();

    updateSystem();
    // TODO: Need to be rewrite to minimize building time accross the solvers.
    // It will be better to directly build all matrices in the solvers.
    makeQPForm();
    sol_->SI_problem(static_cast<int>(c_.size()), constraints_.nrEqConstr, constraints_.nrIneqConstr);

    auto sTime = high_resolution_clock::now();
    bool success = sol_->SI_solve(Q_, c_, Aeq_, beq_, Aineq_, bineq_, lb_, ub_);
    solveTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sTime);

    checkDeleteCostsAndConstraints();
    if (success) {
        updateResults();
    }

    solveAndBuildTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sabTime);
    return success;
}

void LMPC::inform() const noexcept
{
    sol_->SI_inform();
}

double LMPC::solveTime() const noexcept
{
    return solveTime_.count();
}

double LMPC::solveAndBuildTime() const noexcept
{
    return solveAndBuildTime_.count();
}

void LMPC::addCost(const std::shared_ptr<CostFunction>& costFun)
{
    costFun->initializeCost(*ps_);
    spCost_.emplace_back(costFun);
}

void LMPC::addConstraint(const std::shared_ptr<Constraint>& constr)
{
    constr->initializeConstraint(*ps_);
    addConstraintByType(constr);
}

void LMPC::clearCosts() noexcept
{
    spCost_.clear();
}

void LMPC::clearConstraints() noexcept
{
    constraints_.clear();
    clearConstraintMatrices();
}

void LMPC::removeCost(const std::shared_ptr<CostFunction>& costFun)
{
    auto rCost = std::find(spCost_.begin(), spCost_.end(), costFun);
    if (rCost != spCost_.end()) {
        spCost_.erase(rCost);
    }
}

void LMPC::removeConstraint(const std::shared_ptr<Constraint>& constr)
{
    auto rConstr = std::find(constraints_.spConstr.begin(), constraints_.spConstr.end(), constr);
    if (rConstr != constraints_.spConstr.end()) {
        constraints_.spConstr.erase(rConstr);
        switch (constr->constraintType()) {
        case ConstraintFlag::Constraint:
        case ConstraintFlag::EqualityConstraint:
            constraints_.spEqConstr.erase(std::find(constraints_.spEqConstr.begin(), constraints_.spEqConstr.end(), constr));
            break;
        case ConstraintFlag::InequalityConstraint:
            constraints_.spIneqConstr.erase(std::find(constraints_.spIneqConstr.begin(), constraints_.spIneqConstr.end(), constr));
            break;
        case ConstraintFlag::BoundConstraint:
            constraints_.spBoundConstr.erase(std::find(constraints_.spBoundConstr.begin(), constraints_.spBoundConstr.end(), constr));
            break;
        }
    }
}

/*
 *  Protected methods
 */

void LMPC::addConstraintByType(const std::shared_ptr<Constraint>& constr)
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

void LMPC::clearConstraintMatrices()
{
    assert(ps_ != nullptr);

    Aineq_.resize(0, ps_->fullUDim);
    Aeq_.resize(0, ps_->fullUDim);
    bineq_.resize(0);
    beq_.resize(0);
    lb_.setConstant(ps_->fullUDim, -std::numeric_limits<double>::max());
    ub_.setConstant(ps_->fullUDim, std::numeric_limits<double>::max());
}

void LMPC::updateQPMatrixSize()
{
    Q_.resize(ps_->fullUDim, ps_->fullUDim);
    c_.resize(ps_->fullUDim);
}

void LMPC::updateConstraintMatrixSize()
{
    Aeq_.resize(constraints_.nrEqConstr, ps_->fullUDim);
    beq_.resize(constraints_.nrEqConstr);
    Aineq_.resize(constraints_.nrIneqConstr, ps_->fullUDim);
    bineq_.resize(constraints_.nrIneqConstr);
}

void LMPC::updateSystem()
{
    // Reset the QP variables
    Q_.setIdentity();
    Q_ *= 1e-6; // Ensure that Q is positive definite if no cost functions are used
    c_.setZero();

    // Update the system
    if (!ps_->isUpdated) {
        ps_->updateSystem();
    }

    // Update the constraints
    constraints_.updateNr();
    updateConstraintMatrixSize();
    for (auto& cstr : constraints_.spConstr) {
        cstr->update(*ps_);
    }

    // Update the costs
    for (auto& sp : spCost_) {
        sp->update(*ps_);
    }
}

void LMPC::makeQPForm()
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

void LMPC::updateResults()
{
    control_ = sol_->SI_result();
    trajectory_ = ps_->Phi * ps_->x0 + ps_->Psi * control_ + ps_->xi;
}

void LMPC::checkDeleteCostsAndConstraints()
{
    auto check = [](auto& sp, bool useWarn = false, int delLimit = 2) {
        for (auto it = sp.begin(); it != sp.end();) {
            if ((*it).use_count() <= delLimit) {
                CONSTRAINT_DELETION_WARN(useWarn, "%s%s%s", "A '", (*it)->name().c_str(),
                    "' has been destroyed.\nIt has been removed from the controller\n");
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
    check(spCost_, true, 1);
}

} // namespace copra
