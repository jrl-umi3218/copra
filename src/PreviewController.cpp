#include "PreviewController.h"

#include <numeric>
#include "Constrains.h"
#include "PreviewSystem.h"

namespace mpc
{

/*************************************************************************************************
 *                                         MPCTypeFull                                           *
 *************************************************************************************************/

MPCTypeFull::MPCTypeFull(const PreviewSystem &ps, SolverFlag sFlag)
    : nrConstr_(0),
      constr_(),
      sol_(solverFactory(sFlag)),
      Q_(ps.fullUDim, ps.fullUDim),
      AInEq_(3 * ps.fullUDim, 3 * ps.fullUDim), // max size is 1 equality (= 2 inequalities) + 1 inequality
      c_(ps.fullUDim),
      bInEq_(3 * ps.fullUDim),
      Wx_(ps.fullXDim),
      Wu_(ps.fullUDim),
      solveTime_(),
      solveAndBuildTime_()
{
}

void MPCTypeFull::selectQPSolver(SolverFlag flag)
{
    sol_ = solverFactory(flag);
}

bool MPCTypeFull::solve(const PreviewSystem &ps)
{
    solveAndBuildTime_.start();
    makeQPForm(ps);
    sol_->SI_problem(ps.fullUDim, 1, nrConstr_);
    solveTime_.start();
    bool success = sol_->SI_solve(Q_, c_, Eigen::MatrixXd::Zero(0, ps.fullUDim),
                                  Eigen::VectorXd::Zero(0), AInEq_, bInEq_,
                                  Eigen::VectorXd::Constant(ps.fullUDim, -std::numeric_limits<double>::infinity()),
                                  Eigen::VectorXd::Constant(ps.fullUDim, std::numeric_limits<double>::infinity()));
    solveTime_.stop();
    solveAndBuildTime_.stop();
    if (!success)
        sol_->SI_inform();

    return success;
}

const Eigen::VectorXd &MPCTypeFull::control() const noexcept
{
    return sol_->SI_result();
}

Eigen::VectorXd MPCTypeFull::trajectory(const PreviewSystem &ps) const noexcept
{
    return ps.Phi * ps.x0 + ps.Psi * control() + ps.xi;
}

boost::timer::cpu_times MPCTypeFull::solveTime() const noexcept
{
    return solveTime_.elapsed();
}

boost::timer::cpu_times MPCTypeFull::solveAndBuildTime() const noexcept
{
    return solveAndBuildTime_.elapsed();
}

void MPCTypeFull::weights(const PreviewSystem &ps, const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu)
{
    assert(Wx.rows() == ps.xDim);
    assert(Wu.rows() == ps.uDim);

    for (auto i = 0; i < ps.nrStep; ++i)
    {
        Wx_.segment(i * ps.xDim, ps.xDim) = Wx;
        Wu_.segment(i * ps.uDim, ps.uDim) = Wu;
    }
}

void MPCTypeFull::addConstrain(const PreviewSystem &ps, Constrain &constr)
{
    constr_.push_back(&constr);
    constr.initializeConstrain(ps);
    nrConstr_ += constr.nrConstr();

    AInEq_.resize(nrConstr_, ps.fullUDim);
    bInEq_.resize(nrConstr_);
}

void MPCTypeFull::resetConstrains() noexcept
{
    nrConstr_ = 0;
    constr_.clear();
}

/* 
 *  Protected methods
 */

void MPCTypeFull::updateSystem(PreviewSystem &ps)
{
    auto xDim = ps.xDim;
    auto uDim = ps.uDim;

    ps.Phi.block(0, 0, xDim, xDim) = ps.A;
    ps.Psi.block(0, 0, xDim, uDim) = ps.B;
    ps.xi.segment(0, xDim) = ps.d;

    for (auto i = 1; i < ps.nrStep; ++i)
    {
        ps.Phi.block(i * xDim, 0, xDim, xDim) = ps.A * ps.Phi.block((i - 1) * xDim, 0, xDim, xDim);
        for (auto j = 0; j < i; ++j)
        {
            ps.Psi.block(i * xDim, j * uDim, xDim, uDim) = ps.A * ps.Psi.block((i - 1) * xDim, j * uDim, xDim, uDim);
        }
        ps.Psi.block(i * xDim, i * uDim, xDim, uDim) = ps.B;
        ps.xi.segment(i * xDim, xDim) = ps.A * ps.xi.segment((i - 1) * xDim, xDim) + ps.d;
    }

    for (auto cstr : constr_)
        cstr->update(ps);
}

void MPCTypeFull::makeQPForm(const PreviewSystem &ps)
{
    Q_ = ps.Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * ps.Psi + Eigen::MatrixXd(Wu_.asDiagonal());
    c_ = ps.Psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (ps.Phi * ps.x0 - ps.xd + ps.xi);

    int nrLines = 0;
    for (auto cstr : constr_)
    {
        AInEq_.block(nrLines, 0, cstr->nrConstr(), ps.fullUDim) = cstr->A();
        bInEq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }
}

/*************************************************************************************************
 *                                         MPCTypeLast                                           *
 *************************************************************************************************/

MPCTypeLast::MPCTypeLast(const PreviewSystem &ps, SolverFlag sFlag)
    : MPCTypeFull::MPCTypeFull(ps, sFlag)
{
    Wx_.resize(ps.xDim);
}

void MPCTypeLast::weights(const PreviewSystem &ps, const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu)
{
    assert(Wx.rows() == ps.xDim);
    assert(Wu.rows() == ps.uDim);

    Wx_ = Wx;
    for (auto i = 0; i < ps.nrStep; ++i)
        Wu_.segment(i * ps.uDim, ps.uDim) = Wu;
}

/* 
 *  Protected methods
 */

void MPCTypeLast::makeQPForm(const PreviewSystem &ps)
{
    auto xDim = ps.xDim;
    const Eigen::MatrixXd &psi = ps.Psi.bottomRows(xDim);
    Q_ = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * psi + Eigen::MatrixXd(Wu_.asDiagonal());
    c_ = psi.transpose() * Eigen::MatrixXd(Wx_.asDiagonal()) * (ps.Phi.bottomRows(xDim) * ps.x0 - ps.xd.tail(xDim) + ps.xi.tail(xDim));

    int nrLines = 0;
    for (auto cstr : constr_)
    {
        AInEq_.block(nrLines, 0, cstr->nrConstr(), ps.fullUDim) = cstr->A();
        bInEq_.segment(nrLines, cstr->nrConstr()) = cstr->b();
        nrLines += cstr->nrConstr();
    }
}

} // namespace pc