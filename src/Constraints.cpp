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
    if (E.rows() != f.rows())
        throw std::runtime_error("E and f should have same number of rows. E has "
            + std::to_string(E.rows()) + " rows and f has " + std::to_string(f.rows()) + " rows.");
}

void TrajectoryConstraint::trajectory(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    if (E_.rows() != E.rows() || E_.cols() != E.cols())
        throw std::runtime_error("Bad dimension for E. It should be an (" + std::to_string(E_.rows()) + "-by-" + std::to_string(E_.cols())
            + ") matrix but you gave an (" + std::to_string(E.rows()) + "-by-" + std::to_string(E.cols()) + ") matrix");
    if (f_.rows() != f.rows())
        throw std::runtime_error("Bad dimension for f. It should be an (" + std::to_string(f_.rows()) + "-by-1 column vector"
            + ") but you gave an (" + std::to_string(f.rows()) + "-by-1 vector");
    E_ = E;
    f_ = f;
}

void TrajectoryConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (E_.cols() == ps.xDim) {
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;

    } else if (E_.cols() == ps.fullXDim) {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(E_.rows());
    } else {
        throw std::runtime_error("E has a bad dimension. It should be an (nrConstr-by-" + std::to_string(ps.xDim) + ") or (nrConstr-by-"
            + std::to_string(ps.fullXDim) + ") matrix and you gave an (nrConstr-by-" + std::to_string(E_.cols()) + ") matrix");
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
    if (E_.rows() != f.rows())
        throw std::runtime_error("E and f should have same number of rows. E has "
            + std::to_string(E.rows()) + " rows and f has " + std::to_string(f.rows()) + " rows.");
}

void ControlConstraint::control(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    if (E_.rows() != E.rows() || E_.cols() != E.cols())
        throw std::runtime_error("Bad dimension for E. It should be an (" + std::to_string(E_.rows()) + "-by-" + std::to_string(E_.cols())
            + ") matrix but you gave an (" + std::to_string(E.rows()) + "-by-" + std::to_string(E.cols()) + ") matrix");
    if (f_.rows() != f.rows())
        throw std::runtime_error("Bad dimension for f. It should be an (" + std::to_string(f_.rows()) + "-by-1 column vector"
            + ") but you gave an (" + std::to_string(f.rows()) + "-by-1) vector");
    E_ = E;
    f_ = f;
}

void ControlConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (E_.cols() == ps.uDim) {
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrStep;
        A_.resize(nrConstr_, ps.fullUDim);
        b_.resize(nrConstr_);
        A_.setZero();
    } else if (E_.cols() == ps.fullUDim) {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(E_.rows());
        A_ = std::move(E_);
        b_ = std::move(f_);
    } else {
        throw std::runtime_error("E has a bad dimension. It should be an (nrConstr-by-" + std::to_string(ps.uDim) + ") or (nrConstr-by-"
            + std::to_string(ps.fullUDim) + ") matrix but you gave an (nrConstr-by-" + std::to_string(E_.cols()) + ") matrix.");
    }
}

void ControlConstraint::update(const PreviewSystem& ps)
{
    if (!fullSizeEntry_) {
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
    if (lower.rows() != upper.rows())
        throw std::runtime_error("lower and upper should have same number of rows. lower has "
            + std::to_string(lower.rows()) + " rows and upper has " + std::to_string(upper.rows()) + " rows.");

    for (auto line = 0; line < lower_.rows(); ++line) {
        if (lower_(line) != -std::numeric_limits<double>::infinity())
            lowerLines_.push_back(line);
        if (upper_(line) != std::numeric_limits<double>::infinity())
            upperLines_.push_back(line);
    }
}

void TrajectoryBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (lower_.rows() != ps.xDim && lower_.rows() != ps.fullXDim)
        throw std::runtime_error("The lower and upper limit should be of dimension (" + std::to_string(ps.xDim) + "-by-1) or ("
            + std::to_string(ps.fullXDim) + "-by-1) but are of size (" + std::to_string(lower_.rows()) + "-by-1).");

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
    if (lower.rows() != upper.rows())
        throw std::runtime_error("lower and upper should have same number of rows. lower has "
            + std::to_string(lower.rows()) + " rows and upper has " + std::to_string(upper.rows()) + " rows.");
}

void ControlBoundConstraint::controlBound(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
{
    if (lower_.rows() != lower.rows())
        throw std::runtime_error("Bad dimension for lower. It should be an (" + std::to_string(lower_.rows())
            + "-by-1) vector but you gave an (" + std::to_string(lower.rows()) + "-by-1) vector");
    if (upper_.rows() != upper.rows())
        throw std::runtime_error("Bad dimension for upper. It should be an (" + std::to_string(upper_.rows())
            + "-by-1 vector) but you gave an (" + std::to_string(upper.rows()) + "-by-1) vector");
    lower_ = lower;
    upper_ = upper;
}

void ControlBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (lower_.rows() == ps.uDim && upper_.rows() == ps.uDim) {
        nrConstr_ = ps.fullUDim;
        lb_.resize(nrConstr_);
        ub_.resize(nrConstr_);
    } else if (lower_.rows() == ps.fullUDim && upper_.rows() == ps.fullUDim) {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(lower_.rows());
        lb_ = std::move(lower_);
        ub_ = std::move(upper_);
    } else {
        throw std::runtime_error("The lower and upper limit should be of dimension (" + std::to_string(ps.uDim) + "-by-1) or ("
            + std::to_string(ps.fullUDim) + "-by-1) but are of size (" + std::to_string(lower_.rows()) + "-by-1).");
    }
}

void ControlBoundConstraint::update(const PreviewSystem& ps)
{
    if (!fullSizeEntry_) {
        for (int i = 0; i < ps.nrStep; ++i) {
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