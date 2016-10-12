#include "previewSystem.h"

namespace pc
{

PreviewSystemData::PreviewSystemData(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                                     const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                                     std::size_t numberOfSteps)
    : nrStep(numberOfSteps),
      xDim(state.cols()),
      uDim(control.cols()),
      fullXDim(xDim * nrStep),
      fullUDim(uDim * nrStep),
      x0(xInit),
      xd(xTraj),
      A(state),
      B(control),
      d(bias),
      Phi(fullXDim, xDim),
      Psi(fullXDim, fullUDim),
      xi(fullXDim)
{
    Phi.setZero();
    Psi.setZero();
    xi.setZero();
}

System::System(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
               const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
               std::size_t numberOfSteps, PCFlag flag = PCFlag::Last)
    : System::System(state, control, Eigen::VectorXD::Zero(state.rows()), xInit, xTraj, numberOfSteps, flag)
{
}

System::System(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
               const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
               std::size_t numberOfSteps, PCFlag flag = PCFlag::Last)
    : flag_(flag),
      psd_(std::make_unique<PreviewSystemData>(state, control, bias, xInit, xTraj, numberOfSteps)),
      constr_(),
      Q_(psd_->fullUDim, psd_->fullUDim),
      AInEq_(3 * psd_->fullUDim, 3 * psd_->fullUDim), // max size is 1 equality (= 2 inequalities) + 1 inequality
      C_(psd_->fullUDim),
      BInEq_(3 * psd_->fullUDim)
{
    assert(state.row() == control.rows());
    assert(state.rows() == bias.rows());
    assert(state.rows() == xInit.rows());
    assert(state.rows() == xTraj.rows());
}

System::previewSystem()
{
    int i = 0;
    int j = 0;
    auto lineX = [this, i](int var = 0) {
        return (i + var) * xDim_;
    } auto lineU = [this, j](int var = 0) {
        return (j + var) * uDim_;
    }

    for (i = 1; i < nbStep_; ++i)
    {
        psd_->Phi.block(lineX(), 0, xDim_, xDim_) = psd_->A * psd_->Phi..block(lineX(-1), 0, xDim_, xDim_);
        for (j = 0; j < nbStep_; ++j)
            psd_->Psi.block(lineX(), lineU(), xDim_, uDim_) = psd_->A * psd_->Psi.block(lineX(-1), lineU(), xDim_, uDim_);

        psd_->Psi.block(lineX(), lineX(), xDim_, uDim_) = psd_->B;
        psd_->Xi.segment(lineX(), xDim_) = psd_->A * psd_->Xi.segment(lineX(-1), xDim_) + psd_->d;
    }

    for (auto cstr : constr_)
        cstr->update(psd_);
}

void System::makeQPForm()
{
    if (flag_ == PCFlag::Last)
    {
        std::size_t xDim = psd_->xDim;
        Eigen::Matrix psi = psd_->Psi.bottomRows(xDim);
        Q_ = psi.transpose() * psi;
        c_ = 2 * psi.transpose() * (psd_->Phi.bottomRows(xDim) * psd_->x0 - psd_->xd + psd_->xi.tail(xDim));
    }
    else
    {
        Q_ = psd_->Psi.transpose() * psd->Psi;
        c_ = 2 * psd_->Psi.transpose() * (psd - _ > Phi * psd_->x0 - psd_->xd + psd_->xi);
    }

    for (auto i = 0; i < constr_.size(); ++i)
    {
        auto cstr = constr_[î];
        AInEq_.block(i * cstr.nrConstr(), 0, cstr.nrConstr(), psd->fullUDim) = cstr.A();
        bInEq_.segment(i * cstr.nrConstr(), cstr.nrConstr()) = cstr.b();
    }
}

void System::addConstrain(const Constrain &constr)
{
    constr_.push_back(&constr);
    contr.initializeConstrain(psd_);

    std::size_t nrConstrLines = 0;
    for (auto cstr : constr_)
        nrConstrLines += cstr->nrConstr();

    AInEq_.resize(nrConstrLines, psd_->fullUDim);
    bInEq_.resize(nrConstrLines);
}

void System::resetConstrains()
{
    constr_.clear();
}

/**
  * Constrain
  */

Constrain::Constrain(const Eigen::MatrixXd &E, const Eigen::VectorXd f)
    nrConstr_(0),
    E_(E),
    f_(f),
    A_(),
    b_()
{
    assert(E_.rows() == f.rows());
}

void TrajectoryConstrain::initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd)
{
    assert(E_.cols() == psd->xDim);

    nrConstr_ = E_.rows() * psd->nrStep;
    A_.resize(nrConstr_, psd->fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryConstrain::update(const std::unique_ptr<PreviewSystemData> &psd)
{
    std::size_t nrLines = E.row();
    for (std::size_t i = 0; i < psd->nrStep; ++i)
    {
        A_.block(i * nrLines, 0, nrLines, psd->fullUDim) = E_ * psd->Psi.block(i * psd->xDim, 0, psd->xDim, psd->fullUDim);
        b_.segment(i * nrLines, nrLines) = f_ - E_ * (psd->Phi.block(i * psd->xDim, 0, psd->xDim, psd->xDim) * psd->x0 + psd->xi.segment(i * psd->xDim, psd->xDim));
    }
}

void ControlConstrain::initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd)
{
    assert(E_.cols() == psd->uDim);

    nrConstr_ = E_.rows() * psd->nrStep;
    A_.resize(nrConstr_, psd->fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
}

void ControlConstrain::update(const std::unique_ptr<PreviewSystemData> &psd)
{
    std::size_t nrLines = E.row();
    for (std::size_t i = 0; i < psd->nrStep; ++i)
    {
        A_.block(i * nrLines, i * psd->uDim, nrLines, psd->uDim) = E_;
        b_.segment(i * nrLines, nrLines) = f_;
    }
}

} // namespace pc