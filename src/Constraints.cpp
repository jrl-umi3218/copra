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
 *                                         Constraint                                            *
 *************************************************************************************************/

Constraint::Constraint(const std::string& name)
    : name_(name)
    , nrConstr_(0)
{
}

/*************************************************************************************************
 *                           Equality or Inequality Constraint                                   *
 *************************************************************************************************/

EqIneqConstraint::EqIneqConstraint(const std::string& constrQualifier, bool isInequalityConstraint)
    : Constraint(constrQualifier + (isInequalityConstraint ? " inequality constraint" : " equality constraint"))
    , A_()
    , b_()
    , isIneq_(isInequalityConstraint)
{
}

/*************************************************************************************************
 *                                  Trajectory Constraint                                        *
 *************************************************************************************************/

TrajectoryConstraint::TrajectoryConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint)
    : EqIneqConstraint("Trajectory", isInequalityConstraint)
    , E_(E)
    , f_(f)
{
    assert(E_.rows() == f.rows());
}

void TrajectoryConstraint::initializeConstraint(const PreviewSystem& ps)
{
    assert(E_.cols() == ps.xDim);

    nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
    b_.setZero();
}

void TrajectoryConstraint::update(const PreviewSystem& ps)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < ps.nrStep; ++i) {
        A_.block(i * nrLines, 0, nrLines, ps.fullUDim).noalias() = E_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim);
        b_.segment(i * nrLines, nrLines).noalias() = f_ - E_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim));
    }
}

ConstraintFlag TrajectoryConstraint::constraintType()
{
    if (isIneq_)
        return ConstraintFlag::InequalityConstraint;
    else
        return ConstraintFlag::EqualityConstraint;
}

/*************************************************************************************************
 *                                    Control Constraint                                         *
 *************************************************************************************************/

ControlConstraint::ControlConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint)
    : EqIneqConstraint("Control", isInequalityConstraint)
    , E_(E)
    , f_(f)
{
    assert(E_.rows() == f.rows());
}

void ControlConstraint::initializeConstraint(const PreviewSystem& ps)
{
    assert(E_.cols() == ps.uDim);

    nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
    b_.setZero();
}

void ControlConstraint::update(const PreviewSystem& ps)
{
    auto nrLines = static_cast<int>(E_.rows());
    for (int i = 0; i < ps.nrStep; ++i) {
        A_.block(i * nrLines, i * ps.uDim, nrLines, ps.uDim) = E_;
        b_.segment(i * nrLines, nrLines) = f_;
    }
}

ConstraintFlag ControlConstraint::constraintType()
{
    if (isIneq_)
        return ConstraintFlag::InequalityConstraint;
    else
        return ConstraintFlag::EqualityConstraint;
}

/*************************************************************************************************
 *                               Trajectory Bound Constraint                                     *
 *************************************************************************************************/

TrajectoryBoundConstraint::TrajectoryBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
    : EqIneqConstraint("Trajectory bound", true)
    , lower_(lower)
    , upper_(upper)
{
    assert(lower_.rows() == upper_.rows());
}

void TrajectoryBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    assert(lower_.rows() == ps.xDim);

    nrConstr_ = 2 * ps.fullXDim;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
    b_.setZero();
}

void TrajectoryBoundConstraint::update(const PreviewSystem& ps)
{
    A_.block(0, 0, ps.fullXDim, ps.fullUDim).noalias() = ps.Psi;
    A_.block(ps.fullXDim, 0, ps.fullXDim, ps.fullUDim).noalias() = ps.Psi;
    for (int i = 0; i < ps.nrStep; ++i) {
        b_.segment(i * ps.xDim, ps.xDim).noalias() = upper_ - ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 - ps.xi.segment(i * ps.xDim, ps.xDim);
        b_.segment(ps.fullXDim + i * ps.xDim, ps.xDim).noalias() = lower_ - ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 - ps.xi.segment(i * ps.xDim, ps.xDim);
    }
}

ConstraintFlag TrajectoryBoundConstraint::constraintType()
{
    return ConstraintFlag::InequalityConstraint;
}

/*************************************************************************************************
 *                                 Control Bound Constraint                                      *
 *************************************************************************************************/

ControlBoundConstraint::ControlBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
    : Constraint("Control bound constraint")
    , lower_(lower)
    , upper_(upper)
    , lb_()
    , ub_()
{
    assert(lower_.rows() == upper_.rows());
}

void ControlBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    assert(lower_.rows() == ps.uDim);

    nrConstr_ = ps.fullXDim;
    lb_.resize(nrConstr_);
    ub_.resize(nrConstr_);
    lb_.setZero();
    ub_.setZero();
}

void ControlBoundConstraint::update(const PreviewSystem& ps)
{
    for (int i = 0; i < ps.nrStep; ++i) {
        ub_.segment(i * ps.uDim, ps.uDim) = upper_;
        lb_.segment(i * ps.uDim, ps.uDim) = lower_;
    }
}

ConstraintFlag ControlBoundConstraint::constraintType()
{
    return ConstraintFlag::BoundConstraint;
}

} // namespace mpc