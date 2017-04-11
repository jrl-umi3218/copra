// This file is part of mpc.

// mpc is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with mpc.  If not, see
// <http://www.gnu.org/licenses/>.

// header
#include "constraints.h"

// mpc
#include "PreviewSystem.h"
#include "debugUtils.h"

// stl
#include <sstream>

namespace mpc {

/*************************************************************************************************
 *                                         Constraint                                            *
 *************************************************************************************************/

Constraint::Constraint(const std::string& name)
    : name_(name)
    , nrConstr_(0)
    , fullSizeEntry_(false)
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
    checkRows("E", "f", E, f);
}

void TrajectoryConstraint::reset(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    checkMat("E", E, E_);
    checkMat("f", f, f_);

    E_ = E;
    f_ = f;
}

void TrajectoryConstraint::initializeConstraint(const PreviewSystem& ps)
{
    checkColsOnPSXDim("E", E_, &ps);
    if (E_.cols() == ps.xDim) {
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrXStep;

    } else {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(E_.rows());
    }

    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryConstraint::update(const PreviewSystem& ps)
{
    if (fullSizeEntry_) {
        A_.noalias() = E_ * ps.Psi;
        b_.noalias() = f_ - E_ * (ps.Phi * ps.x0 + ps.xi);
    } else {
        auto nrLines = static_cast<int>(E_.rows());
        for (int i = 0; i < ps.nrXStep; ++i) {
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
    checkRows("E", "f", E, f);
}

void ControlConstraint::reset(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    checkMat("E", E, E_);
    checkMat("f", f, f_);

    E_ = E;
    f_ = f;
}

void ControlConstraint::initializeConstraint(const PreviewSystem& ps)
{
    checkColsOnPSUDim("E", E_, &ps);
    if (E_.cols() == ps.uDim) {
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrUStep;
        A_.resize(nrConstr_, ps.fullUDim);
        b_.resize(nrConstr_);
        A_.setZero();
    } else {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(E_.rows());
        A_ = std::move(E_);
        b_ = std::move(f_);
    }
}

void ControlConstraint::update(const PreviewSystem& ps)
{
    if (!fullSizeEntry_) {
        auto nrLines = static_cast<int>(E_.rows());
        for (int i = 0; i < ps.nrUStep; ++i) {
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
 *                                      Mixed Constraint                                         *
 *************************************************************************************************/

MixedConstraint::MixedConstraint(const Eigen::MatrixXd& E, const Eigen::MatrixXd& G, const Eigen::VectorXd& f, bool isInequalityConstraint)
    : EqIneqConstraint("Control", isInequalityConstraint)
    , E_(E)
    , G_(G)
    , f_(f)
{
    checkRows("E", "f", E, f);
    checkRows("G", "f", G, f);
}

void MixedConstraint::reset(const Eigen::MatrixXd& E, const Eigen::MatrixXd& G, const Eigen::VectorXd& f)
{
    checkMat("E", E, E_);
    checkMat("G", G, G_);
    checkMat("f", f, f_);

    E_ = E;
    G_ = G;
    f_ = f;
}

void MixedConstraint::initializeConstraint(const PreviewSystem& ps)
{
    checkColsOnPSXUDim("E", "G", E_, G_, &ps);
    if (E_.cols() == ps.xDim) {
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrXStep;
        A_.resize(nrConstr_, ps.fullUDim);
        b_.resize(nrConstr_);
        A_.setZero();
    } else {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(E_.rows());
        A_ = std::move(E_);
        b_ = std::move(f_);
    }
}

void MixedConstraint::update(const PreviewSystem& /*E*/)
{
    /*
    if (!fullSizeEntry_) {
        auto nrLines = static_cast<int>(E_.rows());
        for (int i = 0; i < ps.nrXStep; ++i) {
            A_.block(i * nrLines, i * ps.uDim, nrLines, ps.uDim) = E_;
            b_.segment(i * nrLines, nrLines) = f_;
        }
    }
    */
}

ConstraintFlag MixedConstraint::constraintType()
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
    checkRows("lower", "upper", lower, upper);

    for (auto line = 0; line < lower_.rows(); ++line) {
        if (lower_(line) != -std::numeric_limits<double>::infinity())
            lowerLines_.push_back(line);
        if (upper_(line) != std::numeric_limits<double>::infinity())
            upperLines_.push_back(line);
    }
}

void TrajectoryBoundConstraint::reset(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
{
    checkMat("lower", lower, lower_);
    checkMat("upper", upper, upper_);

    lower_ = lower;
    upper_ = upper;

    lowerLines_.clear();
    upperLines_.clear();

    for (auto line = 0; line < lower_.rows(); ++line) {
        if (lower_(line) != -std::numeric_limits<double>::infinity())
            lowerLines_.push_back(line);
        if (upper_(line) != std::numeric_limits<double>::infinity())
            upperLines_.push_back(line);
    }
}

void TrajectoryBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    checkRowsOnPSXDim("lower", lower_, &ps);
    if (lower_.rows() == ps.xDim) {
        nrConstr_ = static_cast<int>((lowerLines_.size() + upperLines_.size())) * ps.nrXStep;
    } else {
        nrConstr_ = static_cast<int>((lowerLines_.size() + upperLines_.size()));
        fullSizeEntry_ = true;
    }

    A_.resize(nrConstr_, ps.fullUDim);
    b_.resize(nrConstr_);
}

void TrajectoryBoundConstraint::update(const PreviewSystem& ps)
{
    int nrLines = 0;
    Eigen::VectorXd delta = ps.Phi * ps.x0 + ps.xi;
    for (auto step = 0; step < ps.nrXStep; ++step) {
        for (auto line : lowerLines_) {
            A_.row(nrLines) = ps.Psi.row(line + ps.xDim * step);
            b_(nrLines) = lower_(line) - delta(line + ps.xDim * step);
            ++nrLines;
        }
        if (fullSizeEntry_)
            break;
    }

    for (auto step = 0; step < ps.nrXStep; ++step) {
        for (auto line : upperLines_) {
            A_.row(nrLines) = ps.Psi.row(line + ps.xDim * step);
            b_(nrLines) = upper_(line) - delta(line + ps.xDim * step);
            ++nrLines;
        }
        if (fullSizeEntry_)
            break;
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
    checkRows("lower", "upper", lower, upper);
}

void ControlBoundConstraint::reset(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
{
    checkMat("lower", lower, lower_);
    checkMat("upper", upper, upper_);
    lower_ = lower;
    upper_ = upper;
}

void ControlBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    checkRowsOnPSUDim("lower", lower_, &ps);
    if (lower_.rows() == ps.uDim) {
        nrConstr_ = ps.fullUDim;
        lb_.resize(nrConstr_);
        ub_.resize(nrConstr_);
    } else {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(lower_.rows());
        lb_ = std::move(lower_);
        ub_ = std::move(upper_);
    }
}

void ControlBoundConstraint::update(const PreviewSystem& ps)
{
    if (!fullSizeEntry_) {
        for (int i = 0; i < ps.nrUStep; ++i) {
            ub_.segment(i * ps.uDim, ps.uDim) = upper_;
            lb_.segment(i * ps.uDim, ps.uDim) = lower_;
        }
    }
}

ConstraintFlag ControlBoundConstraint::constraintType()
{
    return ConstraintFlag::BoundConstraint;
}

} // namespace mpc