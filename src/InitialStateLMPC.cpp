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
    : LMPC(sFlag), 
    R_(), 
    r_(),
    E_(), 
    f_(),
    l_(),
    u_(),
    Yeq_(),
    zeq_(),
    Yineq_(),
    zineq_(),
    newQ_(), 
    newc_(), 
    newlb_(), 
    newub_(),
    newAeq_(),
    newbeq_(),
    newAineq_(),
    newbineq_(),
    control_()
{
}

InitialStateLMPC::InitialStateLMPC(const std::shared_ptr<PreviewSystem>& ps, SolverFlag sFlag)
    : LMPC(ps, sFlag), 
    R_(ps_->xDim,ps_->xDim), 
    r_(ps_->xDim),
    E_(ps_->xDim, ps_->fullUDim), 
    f_(ps_->fullUDim),
    l_(ps_->xDim),
    u_(ps_->xDim),
    Yeq_(),
    zeq_(),
    Yineq_(),
    zineq_(),
    newQ_(ps_->xDim+ps_->fullUDim, ps_->xDim+ps_->fullUDim), 
    newc_(ps_->xDim+ps_->fullUDim), 
    newlb_(ps_->xDim+ps_->fullUDim), 
    newub_(ps_->xDim+ps_->fullUDim),
    newAeq_(),
    newbeq_(),
    newAineq_(),
    newbineq_(),
    control_(ps_->fullUDim)
{
    l_.setConstant(-std::numeric_limits<double>::max());
    u_.setConstant(std::numeric_limits<double>::max());
    newlb_.setConstant(-std::numeric_limits<double>::max());
    newub_.setConstant(std::numeric_limits<double>::max());
}


void InitialStateLMPC::initializeController(const std::shared_ptr<PreviewSystem>& ps)
{
    LMPC::initializeController(ps);

    R_.resize(ps_->xDim,ps_->xDim);
    r_.resize(ps_->xDim);
    E_.resize(ps_->xDim, ps_->fullUDim);
    f_.resize(ps_->fullUDim);
    l_.resize(ps_->xDim);
    u_.resize(ps_->xDim);
    newQ_.resize(ps_->xDim+ps_->fullUDim, ps_->xDim+ps_->fullUDim);
    newc_.resize(ps_->xDim+ps_->fullUDim);
    newlb_.resize(ps_->xDim+ps_->fullUDim);
    newub_.resize(ps_->xDim+ps_->fullUDim);
    control_.resize(ps_->fullUDim);

    Eigen::VectorXd xInit = ps->xInit();
    resetInitialStateBounds(xInit,xInit); //this basically prevents optimization of the initial state
    Eigen::MatrixXd R = (1e-6)*Eigen::MatrixXd::Identity(ps_->xDim,ps_->xDim); // ensure that R is positive definite
    Eigen::VectorXd r = Eigen::VectorXd::Zero(ps_->xDim);
    resetInitialStateCost(R, r);
}

bool InitialStateLMPC::solve()
{
    //don't call the function solve() of the base-class!

    using namespace std::chrono;
    auto sabTime = high_resolution_clock::now();

    updateSystem();
    // TODO: Need to be rewrite to minimize building time accross the solvers.
    // It will be better to directly build all matrices in the solvers.
    makeQPForm();
    sol_->SI_problem(ps_->xDim+ps_->fullUDim, constraints_.nrEqConstr, constraints_.nrIneqConstr);

    auto sTime = high_resolution_clock::now();
    bool success = sol_->SI_solve(newQ_, newc_, newAeq_, newbeq_, newAineq_, newbineq_, newlb_, newub_);
    solveTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sTime);

    checkDeleteCostsAndConstraints();

    solveAndBuildTime_ = duration_cast<duration<double>>(high_resolution_clock::now() - sabTime);
    control_ = sol_->SI_result().tail(ps_->fullUDim); //needed here, because of constant function control()
    return success;
}

Eigen::VectorXd InitialStateLMPC::initialState() const noexcept
{
    return sol_->SI_result().head(ps_->xDim);
}

const Eigen::VectorXd& InitialStateLMPC::control() const noexcept
{
    return control_;
}

Eigen::VectorXd InitialStateLMPC::trajectory() const noexcept
{
    return ps_->Phi * initialState() + ps_->Psi * control() + ps_->xi;
}

void InitialStateLMPC::resetInitialStateCost(const Eigen::MatrixXd& R, const Eigen::VectorXd& r)
{
    assert(R_.cols() == R.cols() && R_.rows() == R.rows() && "Matrix size mismatch");
    //R must be positive definite!
    R_ = R;
    r_ = r;
}

void InitialStateLMPC::resetInitialStateBounds(const Eigen::VectorXd& l, const Eigen::VectorXd& u)
{
    assert(l.rows() == ps_->xDim && u.rows() == ps_->xDim && "Vector size mismatch");
    l_ = l;
    u_ = u;
}

/*
 *  Protected methods
 */

void InitialStateLMPC::clearConstraintMatrices()
{
    LMPC::clearConstraintMatrices();

    l_.resize(ps_->xDim);
    u_.resize(ps_->xDim);
    l_.setConstant(-std::numeric_limits<double>::max());
    u_.setConstant(std::numeric_limits<double>::max());
    Yeq_.resize(0, ps_->xDim);
    zeq_.resize(0);
    Yineq_.resize(0, ps_->xDim);
    zineq_.resize(0);
    newAineq_.resize(0, ps_->xDim+ps_->fullUDim);
    newAeq_.resize(0, ps_->xDim+ps_->fullUDim);
    newbineq_.resize(0);
    newbeq_.resize(0);
    newlb_.resize(ps_->xDim+ps_->fullUDim);
    newub_.resize(ps_->xDim+ps_->fullUDim);
    newlb_.setConstant(-std::numeric_limits<double>::max());
    newub_.setConstant(std::numeric_limits<double>::max());
}

void InitialStateLMPC::updateSystem()
{
    LMPC::updateSystem();

    E_.setZero();
    f_.setZero();
    newQ_.setIdentity();
    newQ_ *= 1e-6; // Ensure that newQ is positive definite if no cost functions are used
    newc_.setZero();
    Yeq_.resize(constraints_.nrEqConstr, ps_->xDim);
    zeq_.resize(constraints_.nrEqConstr);
    Yineq_.resize(constraints_.nrIneqConstr, ps_->xDim);
    zineq_.resize(constraints_.nrIneqConstr);
    newAeq_.resize(constraints_.nrEqConstr, ps_->xDim+ps_->fullUDim);
    newbeq_.resize(constraints_.nrEqConstr);
    newAineq_.resize(constraints_.nrIneqConstr, ps_->xDim+ps_->fullUDim);
    newbineq_.resize(constraints_.nrIneqConstr);
}

void InitialStateLMPC::makeQPForm()
{
    //don't call the function makeQPForm() of the base-class!
    
    // Get Costs
    for (auto& cost : spCost_) {
        Q_ += cost->Q();
        c_ += cost->c();
        E_ += cost->E();
        f_ += cost->f();
    }

    int nrLines;
    // Get Equality constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spEqConstr) {
        Aeq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        beq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        Yeq_.block(nrLines, 0, cstr->nrConstr(), ps_->xDim) = cstr->Y();
        zeq_.segment(nrLines, cstr->nrConstr()) = cstr->z();
        nrLines += cstr->nrConstr();
    }
    // Get Inequality constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spIneqConstr) {
        Aineq_.block(nrLines, 0, cstr->nrConstr(), ps_->fullUDim) = cstr->A();
        bineq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        Yineq_.block(nrLines, 0, cstr->nrConstr(), ps_->xDim) = cstr->Y();
        zineq_.segment(nrLines, cstr->nrConstr()) = cstr->z();
        nrLines += cstr->nrConstr();
    }
    // Get Bound constraints
    nrLines = 0;
    for (auto& cstr : constraints_.spBoundConstr) {
        lb_.segment(nrLines, cstr->nrConstr()) = cstr->lower();
        ub_.segment(nrLines, cstr->nrConstr()) = cstr->upper();
        nrLines += cstr->nrConstr();
    }

    // prepare the new QP format
    newQ_.topLeftCorner(ps_->xDim,ps_->xDim) = R_ + E_ * Q_.inverse() * E_.transpose(); //TODO we may want to use another method to compute the inverse of Q
    newQ_.topRightCorner(ps_->xDim,ps_->fullUDim) = E_;
    newQ_.bottomLeftCorner(ps_->fullUDim,ps_->xDim) = E_.transpose();
    newQ_.bottomRightCorner(ps_->fullUDim,ps_->fullUDim) = Q_;
    newc_.head(ps_->xDim) = r_;
    newc_.tail(ps_->fullUDim) = f_;

    newAeq_.leftCols(ps_->xDim) = Yeq_;
    newAeq_.rightCols(ps_->fullUDim) = Aeq_;
    newbeq_ = zeq_;

    newAineq_.leftCols(ps_->xDim) = Yineq_;
    newAineq_.rightCols(ps_->fullUDim) = Aineq_;
    newbineq_ = zineq_;

    newlb_.head(ps_->xDim) = l_;
    newub_.head(ps_->xDim) = u_;
    newlb_.tail(ps_->fullUDim) = lb_;
    newub_.tail(ps_->fullUDim) = ub_;
}

} // namespace copra
