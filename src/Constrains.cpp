#include "Constrains.h"

#include "PreviewSystem.h"

namespace mpc
{

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

void TrajectoryConstrain::initializeConstrain(const PreviewSystem &ps)
{
    assert(E_.cols() == ps.xDim);

    nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryConstrain::update(const PreviewSystem &ps)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < ps.nrStep; ++i)
    {
        A_.block(i * nrLines, 0, nrLines, ps.fullUDim) = E_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim);
        b_.segment(i * nrLines, nrLines) = f_ - E_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim));
    }
}

/*************************************************************************************************
 *                                    Constrol Constrains                                        *
 *************************************************************************************************/

void ControlConstrain::initializeConstrain(const PreviewSystem &ps)
{
    assert(E_.cols() == ps.uDim);

    nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
}

void ControlConstrain::update(const PreviewSystem &ps)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < ps.nrStep; ++i)
    {
        A_.block(i * nrLines, i * ps.uDim, nrLines, ps.uDim) = E_;
        b_.segment(i * nrLines, nrLines) = f_;
    }
}

} // namespace mpc