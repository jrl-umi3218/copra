// This file is part of ModelPreviewController.

// ModelPreviewController is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ModelPreviewController is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ModelPreviewController.  If not, see
// <http://www.gnu.org/licenses/>.

#include "Constraints.h"

#include "PreviewSystem.h"

namespace mpc {

/*************************************************************************************************
 *                                          Constrains                                           *
 *************************************************************************************************/

Constraint::Constraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
    : nrConstr_(0)
    , E_(E)
    , A_()
    , f_(f)
    , b_()
{
    assert(E_.rows() == f.rows());
}

/*************************************************************************************************
 *                                  Trajectory Constrains                                        *
 *************************************************************************************************/

void TrajectoryConstraint::initializeConstraint(const PreviewSystem& ps)
{
    assert(E_.cols() == ps.xDim);

    nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryConstraint::update(const PreviewSystem& ps)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < ps.nrStep; ++i) {
        A_.block(i * nrLines, 0, nrLines, ps.fullUDim) = E_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim);
        b_.segment(i * nrLines, nrLines) = f_ - E_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim));
    }
}

std::string TrajectoryConstraint::name() const noexcept
{
    return "Trajectory constraint";
}

/*************************************************************************************************
 *                                    Constrol Constrains                                        *
 *************************************************************************************************/

void ControlConstraint::initializeConstraint(const PreviewSystem& ps)
{
    assert(E_.cols() == ps.uDim);

    nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
}

void ControlConstraint::update(const PreviewSystem& ps)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < ps.nrStep; ++i) {
        A_.block(i * nrLines, i * ps.uDim, nrLines, ps.uDim) = E_;
        b_.segment(i * nrLines, nrLines) = f_;
    }
}

std::string ControlConstraint::name() const noexcept
{
    return "Control constraint";
}

} // namespace mpc