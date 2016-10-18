#include "previewController.h"

#include <numeric>

namespace pc
{

/*************************************************************************************************
 *                                     Preview System Data                                       *
 *************************************************************************************************/

PreviewSystemData::PreviewSystemData(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                                     const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                                     int numberOfSteps)
    : nrStep(numberOfSteps),
      xDim(static_cast<int>(state.cols())),
      uDim(static_cast<int>(control.cols())),
      fullXDim(xDim * nrStep),
      fullUDim(uDim * nrStep), x0(xInit),
      xd(fullXDim),
      A(state),
      B(control),
      d(bias),
      Phi(fullXDim, xDim),
      Psi(fullXDim, fullUDim),
      xi(fullXDim)
{
    assert(xTraj.rows() == xDim || xTraj.rows() == fullXDim);

    auto xTrajLen = static_cast<int>(xTraj.rows());
    for (auto i = 0; i < fullXDim; i += xTrajLen)
        xd.segment(i, xTrajLen) = xTraj;
    Phi.setZero();
    Psi.setZero();
    xi.setZero();
}

/*************************************************************************************************
 *                                     Preview Controller                                        *
 *************************************************************************************************/

PreviewController::PreviewController(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                                     const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                                     int numberOfSteps, PCFlag pcFlag, SolverFlag sFlag)
    : PreviewController::PreviewController(state, control, Eigen::VectorXd::Zero(state.rows()), xInit, xTraj, numberOfSteps, pcFlag, sFlag)
{
}

PreviewController::PreviewController(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                                     const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                                     int numberOfSteps, PCFlag pcFlag, SolverFlag sFlag)
    : pcFlag_(pcFlag),
      nrConstr_(0),
      psd_(std::make_unique<PreviewSystemData>(state, control, bias, xInit, xTraj, numberOfSteps)),
      constr_(),
      sol_(solverFactory(sFlag)),
      Q_(psd_->fullUDim, psd_->fullUDim),
      AInEq_(3 * psd_->fullUDim, 3 * psd_->fullUDim), // max size is 1 equality (= 2 inequalities) + 1 inequality
      c_(psd_->fullUDim),
      bInEq_(3 * psd_->fullUDim),
      Wx_(pcFlag == PCFlag::Last ? psd_->xDim : psd_->fullXDim),
      Wu_(psd_->fullUDim),
      solveTime_(),
      solveAndBuildTime_()
{
    assert(state.rows() == control.rows());
    assert(state.rows() == bias.rows());
    assert(state.rows() == xInit.rows());
    assert(state.rows() == xTraj.rows());
    assert(state.cols() == xInit.rows());
}

void PreviewController::selectQPSolver(SolverFlag flag)
{
    sol_ = solverFactory(flag);
}

bool PreviewController::solve()
{
    solveAndBuildTime_.start();
    previewSystem();
    makeQPForm();
    sol_->SI_problem(psd_->fullUDim, 1, nrConstr_);
    solveTime_.start();
    bool success = sol_->SI_solve(Q_, c_, Eigen::MatrixXd::Zero(0, psd_->fullUDim),
                                  Eigen::VectorXd::Zero(0), AInEq_, bInEq_,
                                  Eigen::VectorXd::Constant(psd_->fullUDim, -std::numeric_limits<double>::infinity()),
                                  Eigen::VectorXd::Constant(psd_->fullUDim, std::numeric_limits<double>::infinity()));
    solveTime_.stop();
    solveAndBuildTime_.stop();
    if(!success)
        sol_->SI_inform();

    return success;
}

const Eigen::VectorXd &PreviewController::control() const noexcept
{
    return sol_->SI_result();
}

Eigen::VectorXd PreviewController::trajectory() const noexcept
{
    return psd_->Phi * psd_->x0 + psd_->Psi * control() + psd_->xi;
}

boost::timer::cpu_times PreviewController::solveTime() const noexcept
{
    return solveTime_.elapsed();
}

boost::timer::cpu_times PreviewController::solveAndBuildTime() const noexcept
{
    return solveAndBuildTime_.elapsed();
}

void PreviewController::weights(double wx, double wu)
{
    Wx_.fill(wx);
    Wu_.fill(wu);
}

void PreviewController::weights(const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu)
{
    assert(Wx.rows() == psd_->xDim);
    assert(Wu.rows() == psd_->uDim);

    if (pcFlag_ == PCFlag::Last)
        Wx_ = Wx;
    else
        for (auto i = 0; i < psd_->nrStep; ++i)
            Wx_.segment(i * psd_->xDim, psd_->xDim) = Wx;

    for (auto i = 0; i < psd_->nrStep; ++i)
        Wu_.segment(i * psd_->uDim, psd_->uDim) = Wu;
}

void PreviewController::addConstrain(Constrain &constr)
{
    constr_.push_back(&constr);
    constr.initializeConstrain(psd_);
    nrConstr_ += constr.nrConstr();

    AInEq_.resize(nrConstr_, psd_->fullUDim);
    bInEq_.resize(nrConstr_);
}

void PreviewController::resetConstrains() noexcept
{
    nrConstr_ = 0;
    constr_.clear();
}

/* 
 *  Private methods
 */

void PreviewController::previewSystem()
{
    auto xDim = psd_->xDim;
    auto uDim = psd_->uDim;

    psd_->Phi.block(0, 0, xDim, xDim) = psd_->A;
    psd_->Psi.block(0, 0, xDim, uDim) = psd_->B;
    psd_->xi.segment(0, xDim) = psd_->d;

    for (auto i = 1; i < psd_->nrStep; ++i)
    {
        psd_->Phi.block(i * xDim, 0, xDim, xDim) = psd_->A * psd_->Phi.block((i - 1) * xDim, 0, xDim, xDim);
        for (auto j = 0; j < i; ++j)
        {
            psd_->Psi.block(i * xDim, j * uDim, xDim, uDim) = psd_->A * psd_->Psi.block((i - 1) * xDim, j * uDim, xDim, uDim);
        }
        psd_->Psi.block(i * xDim, i * uDim, xDim, uDim) = psd_->B;
        psd_->xi.segment(i * xDim, xDim) = psd_->A * psd_->xi.segment((i - 1) * xDim, xDim) + psd_->d;
    }

    for (auto cstr : constr_)
        cstr->update(psd_);
}

void PreviewController::makeQPForm()
{
    if (pcFlag_ == PCFlag::Last)
    {
        auto xDim = psd_->xDim;
        const Eigen::MatrixXd &psi = psd_->Psi.bottomRows(xDim);
        Q_ = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * psi + Eigen::MatrixXd(Wu_.asDiagonal());
        c_ = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (psd_->Phi.bottomRows(xDim) * psd_->x0 - psd_->xd.tail(xDim) + psd_->xi.tail(xDim));
    }
    else
    {
        Q_ = psd_->Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * psd_->Psi + Eigen::MatrixXd(Wu_.asDiagonal());
        c_ = psd_->Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (psd_->Phi * psd_->x0 - psd_->xd + psd_->xi);
    }
    int nrLines = 0;
    for (auto cstr : constr_)
    {
        AInEq_.block(nrLines, 0, cstr->nrConstr(), psd_->fullUDim) = cstr->A();
        bInEq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }
}

/*************************************************************************************************
 *                                          Constrains                                           *
 *************************************************************************************************/

Constrain::Constrain(const Eigen::MatrixXd &E, const Eigen::VectorXd &f)
    : nrConstr_(0),
      E_(E),
      A_(),
      f_(f),
      b_()
{
    assert(E_.rows() == f.rows());
}

/*************************************************************************************************
 *                                  Trajectory Constrains                                        *
 *************************************************************************************************/

void TrajectoryConstrain::initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd)
{
    assert(E_.cols() == psd->xDim);

    nrConstr_ = static_cast<int>(E_.rows()) * psd->nrStep;
    A_.resize(nrConstr_, psd->fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryConstrain::update(const std::unique_ptr<PreviewSystemData> &psd)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < psd->nrStep; ++i)
    {
        A_.block(i * nrLines, 0, nrLines, psd->fullUDim) = E_ * psd->Psi.block(i * psd->xDim, 0, psd->xDim, psd->fullUDim);
        b_.segment(i * nrLines, nrLines) = f_ - E_ * (psd->Phi.block(i * psd->xDim, 0, psd->xDim, psd->xDim) * psd->x0 + psd->xi.segment(i * psd->xDim, psd->xDim));
    }
}

/*************************************************************************************************
 *                                    Constrol Constrains                                        *
 *************************************************************************************************/

void ControlConstrain::initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd)
{
    assert(E_.cols() == psd->uDim);

    nrConstr_ = static_cast<int>(E_.rows()) * psd->nrStep;
    A_.resize(nrConstr_, psd->fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
}

void ControlConstrain::update(const std::unique_ptr<PreviewSystemData> &psd)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < psd->nrStep; ++i)
    {
        A_.block(i * nrLines, i * psd->uDim, nrLines, psd->uDim) = E_;
        b_.segment(i * nrLines, nrLines) = f_;
    }
}

} // namespace pc