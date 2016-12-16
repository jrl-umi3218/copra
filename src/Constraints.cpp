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
#include <iostream>

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
    assert(E_.cols() == ps.uDim || E_.cols() == ps.fullUDim);

    if (E_.cols() == ps.uDim)
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    else
        nrConstr_ = static_cast<int>(E_.rows());

    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryConstraint::update(const PreviewSystem& ps)
{
    if (E_.rows() == nrConstr_) {
        A_ = E_ * ps.Psi;
        b_ = f_ - E_ * (ps.Phi * ps.x0 + ps.xi);
    } else {
        auto nrLines = static_cast<int>(E_.rows());
        for (int i = 0; i < ps.nrStep; ++i) {
            A_.block(i * nrLines, 0, nrLines, ps.fullUDim).noalias() = E_ * ps.Psi.block(i * ps.xDim, 0, ps.xDim, ps.fullUDim);
            b_.segment(i * nrLines, nrLines).noalias() = f_ - E_ * (ps.Phi.block(i * ps.xDim, 0, ps.xDim, ps.xDim) * ps.x0 + ps.xi.segment(i * ps.xDim, ps.xDim));
        }
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
    std::cout << E_.cols() << " " << ps.fullUDim << std::endl;
    assert(E_.cols() == ps.uDim || E_.cols() == ps.fullUDim);

    if (E_.cols() == ps.uDim)
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
    else
        nrConstr_ = static_cast<int>(E_.rows());

    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
    A_.setZero();
}

void ControlConstraint::update(const PreviewSystem& ps)
{
    if (E_.rows() == nrConstr_) {
        A_ = E_;
        b_ = f_;
    } else {
        auto nrLines = static_cast<int>(E_.rows());
        for (int i = 0; i < ps.nrStep; ++i) {
            A_.block(i * nrLines, i * ps.uDim, nrLines, ps.uDim) = E_;
            b_.segment(i * nrLines, nrLines) = f_;
        }
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
    , lowerLines_()
    , upperLines_()
{
    assert(lower_.rows() == upper_.rows());
    for (auto line = 0; line < lower_.rows(); ++line) {
        if (lower_(line) != -std::numeric_limits<double>::infinity())
            lowerLines_.push_back(line);
        if (upper_(line) != std::numeric_limits<double>::infinity())
            upperLines_.push_back(line);
    }
}

void TrajectoryBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    assert(lower_.rows() == ps.xDim);

    nrConstr_ = static_cast<int>((lowerLines_.size() + upperLines_.size())) * ps.nrStep;
    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryBoundConstraint::update(const PreviewSystem& ps)
{
    int nrLines = 0;
    Eigen::VectorXd delta = ps.Phi * ps.x0 + ps.xi;
    for (auto step = 0; step < ps.nrStep; ++step) {
        for (auto line : lowerLines_) {
            A_.row(nrLines) = ps.Psi.row(line + ps.xDim * step);
            b_(nrLines) = lower_(line) - delta(line + ps.xDim * step);
            ++nrLines;
        }
    }

    for (auto step = 0; step < ps.nrStep; ++step) {
        for (auto line : upperLines_) {
            A_.row(nrLines) = ps.Psi.row(line + ps.xDim * step);
            b_(nrLines) = upper_(line) - delta(line + ps.xDim * step);
            ++nrLines;
        }
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

    nrConstr_ = ps.fullUDim;
    lb_.resize(nrConstr_);
    ub_.resize(nrConstr_);
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