/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "InitialStateLMPC.h"
#include "constraints.h"
#include "costFunctions.h"

namespace copra {

/*************************************************************************************************
 *                                             InitialStateLMPC                                               *
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

// // TODO: Better split base class function (the SI_Solve) to not have redundant code
// bool InitialStateLMPC::solve()
// {
//     using namespace std::chrono;
//     auto sabTime = high_resolution_clock::now();

//     updateSystem();
//     // TODO: Need to be rewrite to minimize building time accross the solvers.
//     // It will be better to directly build all matrices in the solvers.
//     makeQPForm();
//     sol_->SI_problem(ps_->xDim + ps_->fullUDim, constraints_.nrEqConstr, constraints_.nrIneqConstr);

//     auto sTime = high_resolution_clock::now();
//     bool success = sol_->SI_solve(newQ_, newc_, newAeq_, newbeq_, newAineq_, newbineq_, newlb_, newub_);
//     solveTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sTime);

//     checkDeleteCostsAndConstraints();

//     solveAndBuildTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sabTime);
//     control_ = sol_->SI_result().tail(ps_->fullUDim); //needed here, because of constant function control()
//     return success;
// }

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
    // Final QP matrix resize
    // newAineq_.resize(0, ps_->xDim + ps_->fullUDim);
    // newAeq_.resize(0, ps_->xDim + ps_->fullUDim);
    // newbineq_.resize(0);
    // newbeq_.resize(0);
    // newlb_.resize(ps_->xDim + ps_->fullUDim);
    // newub_.resize(ps_->xDim + ps_->fullUDim);
    // newlb_.setConstant(-std::numeric_limits<double>::max());
    // newub_.setConstant(std::numeric_limits<double>::max());
    int optSize = ps_->xDim + ps_->fullUDim;
    Aineq_.resize(0, optSize);
    Aeq_.resize(0, optSize);
    bineq_.resize(0);
    beq_.resize(0);
    lb_.setConstant(optSize, -std::numeric_limits<double>::max());
    ub_.setConstant(optSize, std::numeric_limits<double>::max());
    // Base QP matrix resize
    // Aineq_.resize(0, ps_->fullUDim);
    // Aeq_.resize(0, ps_->fullUDim);
    // bineq_.resize(0);
    // beq_.resize(0);
    // lb_.resize(ps_->fullUDim);
    // ub_.resize(ps_->fullUDim);
    // lb_.setConstant(-std::numeric_limits<double>::max());
    // ub_.setConstant(std::numeric_limits<double>::max());
    // l_.setConstant(ps_->xDim, -std::numeric_limits<double>::max());
    // u_.setConstant(ps_->xDim, std::numeric_limits<double>::max());
    // Yeq_.resize(0, ps_->xDim);
    // zeq_.resize(0);
    // Yineq_.resize(0, ps_->xDim);
    // zineq_.resize(0);
    // baseAineq_.resize(0, ps_->fullUDim);
    // baseAeq_.resize(0, ps_->fullUDim);
    // basebineq_.resize(0);
    // basebeq_.resize(0);
    // baselb_.setConstant(ps_->fullUDim, -std::numeric_limits<double>::max());
    // baseub_.setConstant(ps_->fullUDim, std::numeric_limits<double>::max());
}

void InitialStateLMPC::updateQPMatrix()
{
    // Final QP matrix resize
    // newQ_.resize(ps_->xDim + ps_->fullUDim, ps_->xDim + ps_->fullUDim);
    // newc_.resize(ps_->xDim + ps_->fullUDim);
    Q_.resize(ps_->xDim + ps_->fullUDim, ps_->xDim + ps_->fullUDim);
    c_.resize(ps_->xDim + ps_->fullUDim);

    // Base QP matrix resize
    // E_.resize(ps_->xDim, ps_->fullUDim);
    // f_.resize(ps_->fullUDim);
    // Q_.resize(ps_->fullUDim, ps_->fullUDim);
    // c_.resize(ps_->fullUDim);
    // baseQ_.resize(ps_->fullUDim, ps_->fullUDim);
    // basec_.resize(ps_->fullUDim);
    // newlb_.resize(ps_->xDim + ps_->fullUDim);
    // newub_.resize(ps_->xDim + ps_->fullUDim);
    // control_.resize(ps_->fullUDim);

    // Eigen::VectorXd xInit = ps->x0;
    // resetInitialStateBounds(ps->x0, ps->x0); //this basically prevents optimization of the initial state
    // Eigen::MatrixXd R = Eigen::MatrixXd::Identity(ps_->xDim, ps_->xDim) * 1e-6; // ensure that R is positive definite
    // Eigen::VectorXd r = Eigen::VectorXd::Zero(ps_->xDim);
    // resetInitialStateCost(R, r);
    // l_ = ps_->x0; // This basically prevents optimization of the initial state
    // u_ = ps_->x0; // This basically prevents optimization of the initial state
    // R_ = Eigen::MatrixXd::Identity(ps_->xDim, ps_->xDim) * 1e-6; // ensure that R is positive definite
    // r_.setZero(ps_ - xDim);
}

void InitialStateLMPC::updateConstraintMatrixSize()
{
    // E_.setZero();
    // f_.setZero();
    // baseQ_.setIdentity();
    // baseQ_ *= 1e-6; // Ensure that newQ is positive definite if no cost functions are used
    // basec_.setZero();
    // Yeq_.resize(constraints_.nrEqConstr, ps_->xDim);
    // zeq_.resize(constraints_.nrEqConstr);
    // Yineq_.resize(constraints_.nrIneqConstr, ps_->xDim);
    // zineq_.resize(constraints_.nrIneqConstr);
    Aeq_.resize(constraints_.nrEqConstr, ps_->xDim + ps_->fullUDim);
    beq_.resize(constraints_.nrEqConstr);
    Aineq_.resize(constraints_.nrIneqConstr, ps_->xDim + ps_->fullUDim);
    bineq_.resize(constraints_.nrIneqConstr);
    // baseAeq_.resize(constraints_.nrEqConstr, ps_->fullUDim);
    // basebeq_.resize(constraints_.nrEqConstr);
    // baseAineq_.resize(constraints_.nrIneqConstr, ps_->fullUDim);
    // basebineq_.resize(constraints_.nrIneqConstr);
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

    // prepare the new QP format
    // newQ_.topLeftCorner(ps_->xDim, ps_->xDim) = R_ + E_ * Q_.inverse() * E_.transpose(); //TODO we may want to use another method to compute the inverse of Q
    // newQ_.topRightCorner(ps_->xDim, ps_->fullUDim) = E_;
    // newQ_.bottomLeftCorner(ps_->fullUDim, ps_->xDim) = E_.transpose();
    // newQ_.bottomRightCorner(ps_->fullUDim, ps_->fullUDim) = Q_;
    // newc_.head(ps_->xDim) = r_;
    // newc_.tail(ps_->fullUDim) = f_;

    // newAeq_.leftCols(ps_->xDim) = Yeq_;
    // newAeq_.rightCols(ps_->fullUDim) = Aeq_;
    // newbeq_ = zeq_;

    // newAineq_.leftCols(ps_->xDim) = Yineq_;
    // newAineq_.rightCols(ps_->fullUDim) = Aineq_;
    // newbineq_ = zineq_;

    // newlb_.head(ps_->xDim) = l_;
    // newub_.head(ps_->xDim) = u_;
    // newlb_.tail(ps_->fullUDim) = lb_;
    // newub_.tail(ps_->fullUDim) = ub_;
}

void InitialStateLMPC::updateResults()
{
    control_ = sol_->SI_result().tail(ps_->fullUDim);
    trajectory_ = ps_->Phi * initialState() + ps_->Psi * control_ + ps_->xi;
}

} // namespace copra
