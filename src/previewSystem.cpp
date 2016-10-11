#include "previewSystem.h"

namespace pc
{

PreviewSystemData::PreviewSystemData(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
        const Eigen::VectorXd& d, std::size_t nbStep) :
    nbStep(nbStep),
    xDim(A.cols()),
    uDim(B.cols()),
    fullXDim(xDim_*nbStep_),
    fullUDim(uDim_*nbStep_),
    A(A),
    B(B),
    d(),
    Phi(fullXDim_, xDim_),
    Psi(fullXDim_, fullUDim_),
    xi(fullXDim_)
{
    Phi.setZero();
    Psi.setZero();
    xi.setZero();
}


System::System(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::VectorXd& d,
        std::size_t nbStep, PCFlag flag) :
    flag_(flag),
    psd_(std::make_unique<PreviewSystem>(A, B, d, nbStep)),
    constr_()
    Q_(fullUDim_, fullUDim_),
    AInEq_(3*fullUDim_, 3*fullUDim_), // max size is 1 equality (= 2 inequalities) + 1 inequality
    C_(fullUDim_),
    BInEq_(3*fullUDim_)
{
    assert(A.row() == B.rows());
    assert(A.rows() == d.rows());
}

System::previewSystem()
{
    int i = 0;
    int j = 0;
    auto lineX = [this, i](int var=0)
    {
        return (i + var)*xDim_;
    }
    auto lineU = [this, j](int var=0)
    {
        return (j + var)*uDim_;
    }

    for(i = 1; i < nbStep_; ++i)
    {
        psd_->Phi.block(lineX(), 0, xDim_, xDim_) = psd_->A*psd_->Phi..block(lineX(-1), 0, xDim_, xDim_);
        for(j = 0; j < nbStep_; ++j)
            psd_->Psi.block(lineX(), lineU(), xDim_, uDim_) = psd_->A*psd_->Psi.block(lineX(-1), lineU(), xDim_, uDim_)

        psd_->Psi.block(lineX(), lineX(), xDim_, uDim_) = psd_->B;
        psd_->Xi.segment(lineX(), xDim_) = psd_->A*psd_->Xi.segment(lineX(-1), xDim_) + psd_->d; 
    }

    for(auto cstr: constr_)
        cstr->update(psd_);
}

void System::makeQPForm()
{
    std::size_t nrConstrLines = 0;
    for(auto cstr: constr_)
        nrConstrLines += constr_->A().rows();

    AInEq_.resize(nrConstrLines, psd_->fullUDim);
    bInEq_.resize(nrConstrLines);

    if(flag_ == PCFlag::Last)
    {
        std::size_t xDim = psd_->xDim;
        Eigen::Matrix psi = psd_->Psi.bottomRows(xDim)
        Q_ = psi.transpose()*psi;
        c_ = 2*psi.transpose()*(psd_->Phi.bottomRows(xDim)*psd_->x0 - psd_->xd + psd_->xi.tail(xDim));
    }
    else
    {
        Q_ = psd_->Psi.transpose()*psd->Psi;
        c_ = 2*psd_->Psi.transpose()*(psd-_>Phi*psd_->x0 - psd_->xd + psd_->xi);
    }

    for(auto i = 0; i < constr_.size(); ++i)
    {
        auto cstr = constr_[î];
        AInEq_.block(i*cstr.nrConstr(), 0, cstr.nrConstr(), psd->fullUDim) = cstr.A();
        bInEq_.segment(i*cstr.nrConstr(), cstr.nrConstr()) = cstr.b();
    }
}

void System::addConstrain(const Constrain& constr)
{
    constr_.push_back(&constr);
    contr.initializeConstrain(psd_);
}

void System::resetConstrains()
{
    constr_.clear();
}


/**
  * Constrain
  */


Constrain::Constrain(const Eigen::MatrixXd& E, const Eigen::VectorXd f)
    nrConstr_(0),
    E_(E),
    f_(f),
    A_(),
    b_()
{
    assert(E_.rows() == f.rows());
}

void TrajectoryConstrain::initializeConstrain(const std::unique_ptr<PreviewSystemData>& psd)
{
    assert(E_.cols() == psd->xDim);

    nrConstr_ = E_.rows()*psd->nbStep;
    A_.resize(nrConstr_, psd->fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryConstrain::update(const std::unique_ptr<PreviewSystemData>& psd)
{
    std::size_t nrLines = E.row();
    for(std::size_t i = 0; i < psd->nbStep; ++i)
    {
        A_.block(i*nrLines, 0, nrLines, psd->fullUDim) = E_*psd->Psi.block(i*psd->xDim, 0, psd->xDim, psd->fullUDim);
        b_.segment(i*nrLines, nrLines) = f_ - E_*(psd->Phi.block(i*psd->xDim, 0, psd->xDim, psd->xDim)*psd->x0 \
            + psd->xi.segment(i*psd->xDim, psd->xDim));
    }
}

void ControlConstrain::initializeConstrain(const std::unique_ptr<PreviewSystemData>& psd)
{
    assert(E_.cols() == psd->uDim);

    nrConstr_ = E_.rows()*psd->nbStep;
    A_.resize(nrConstr_, psd->fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
}

void ControlConstrain::update(const std::unique_ptr<PreviewSystemData>& psd)
{
    std::size_t nrLines = E.row();
    for(std::size_t i = 0; i < psd->nbStep; ++i)
    {
        A_.block(i*nrLines, i*psd->uDim, nrLines, psd->uDim) = E_;
        b_.segment(i*nrLines, nrLines) = f_;
    }
}

} // namespace pc