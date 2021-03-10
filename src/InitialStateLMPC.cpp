/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "InitialStateLMPC.h"
#include "constraints.h"
#include "costFunctions.h"

namespace copra {

/*************************************************************************************************
 *                                             InitialStateLMPC                                  *
 *************************************************************************************************/

InitialStateLMPC::InitialStateLMPC(SolverFlag sFlag)
    : LMPC(sFlag)
{
}

InitialStateLMPC::InitialStateLMPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag)
    : LMPC(sFlag)
    , R_(Eigen::MatrixXd::Zero(ps->xDim, ps->xDim))
    , r_(Eigen::VectorXd::Zero(ps->xDim))
    , x0lb_(ps->x0)
    , x0ub_(ps->x0)
{
    initializeController(ps);
}

Eigen::VectorXd InitialStateLMPC::initialState() const noexcept
{
    return sol_->SI_result().head(ps_->xDim);
}

void InitialStateLMPC::resetInitialStateCost(const Eigen::MatrixXd& R, const Eigen::VectorXd& r)
{
    // R must be positive definite!
    R_ = R;
    r_ = r;
}

void InitialStateLMPC::resetInitialStateBounds(const Eigen::VectorXd& x0lb, const Eigen::VectorXd& x0ub)
{
    x0lb_ = x0lb;
    x0ub_ = x0ub;
}

/*
 *  Protected methods
 */

void InitialStateLMPC::clearConstraintMatrices()
{
    int optSize = ps_->xDim + ps_->fullUDim;
    Aineq_.resize(0, optSize);
    Aeq_.resize(0, optSize);
    bineq_.resize(0);
    beq_.resize(0);
    lb_.setConstant(optSize, -std::numeric_limits<double>::max());
    ub_.setConstant(optSize, std::numeric_limits<double>::max());
}

void InitialStateLMPC::updateQPMatrixSize()
{
    Q_.resize(ps_->xDim + ps_->fullUDim, ps_->xDim + ps_->fullUDim);
    c_.resize(ps_->xDim + ps_->fullUDim);

void InitialStateLMPC::updateConstraintMatrixSize()
{
    Aeq_.resize(constraints_.nrEqConstr, ps_->xDim + ps_->fullUDim);
    beq_.resize(constraints_.nrEqConstr);
    Aineq_.resize(constraints_.nrIneqConstr, ps_->xDim + ps_->fullUDim);
    bineq_.resize(constraints_.nrIneqConstr);
}

void InitialStateLMPC::makeQPForm()
{
    // Get Costs
    for (auto& cost : spCost_) {
        Q_.bottomRightCorner(ps_->fullUDim, ps_->fullUDim) += cost->Q();
        Q_.topRightCorner(ps_->xDim, ps_->fullUDim) += cost->E();
        c_.tail(ps_->fullUDim) += cost->f();
    }

    int nrLines = 0;
    // Get Equality constraints
    for (auto& cstr : constraints_.spEqConstr) {
        Aeq_.block(nrLines, 0, cstr->nrConstr(), ps_->xDim) = cstr->Y();
        Aeq_.block(nrLines, ps_->xDim, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        beq_.segment(nrLines, cstr->nrConstr()) = cstr->z();
        nrLines += cstr->nrConstr();
    }

    // Get Inequality constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spIneqConstr) {
        Aineq_.block(nrLines, 0, cstr->nrConstr(), ps_->xDim) = cstr->Y();
        Aineq_.block(nrLines, ps_->xDim, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        bineq_.segment(nrLines, cstr->nrConstr()) = cstr->z();
        nrLines += cstr->nrConstr();
    }

    // Get Bound constraints
    nrLines = ps_->xDim;
    for (auto& cstr : constraints_.spBoundConstr) {
        lb_.segment(nrLines, cstr->nrConstr()) = cstr->lower();
        ub_.segment(nrLines, cstr->nrConstr()) = cstr->upper();
        nrLines += cstr->nrConstr();
    }

    // Update last part of Q
    {
        auto E = Q_.topRightCorner(ps_->xDim, ps_->fullUDim); // Carefull! Temporary Eigen Proxy class
        auto Q = Q_.bottomRightCorner(ps_->fullUDim, ps_->fullUDim); // Carefull! Temporary Eigen Proxy class
        Q_.bottomLeftCorner(ps_->fullUDim, ps_->xDim) = E.transpose();
        Q_.topLeftCorner(ps_->xDim, ps_->xDim).noalias() = R_ + E * Q.inverse() * E.transpose(); //TODO we may want to use another method to compute the inverse of Q
    }
    c_.head(ps_->xDim) = r_;
    lb_.head(ps_->xDim) = x0lb_;
    ub_.head(ps_->xDim) = x0ub_;
}

void InitialStateLMPC::updateResults()
{
    control_ = sol_->SI_result().tail(ps_->fullUDim);
    trajectory_ = ps_->Phi * initialState() + ps_->Psi * control_ + ps_->xi;
}

} // namespace copra
