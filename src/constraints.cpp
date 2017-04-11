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
    if (E.rows() != f.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("E", "f", E, f));
}

void TrajectoryConstraint::reset(const Eigen::MatrixXd& E, const Eigen::VectorXd& f)
{
    if (E.rows() != E_.rows() || E.cols() != E_.cols())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("E", E, E_));
    if (f.rows() != f_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("f", f, f_));

    E_ = E;
    f_ = f;
}

void TrajectoryConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (E_.cols() == ps.xDim) {
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrXStep;

    } else if (E_.cols() == ps.fullXDim) {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(E_.rows());
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXDim("E", E_, &ps));
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

ControlConstraint::ControlConstraint(const Eigen::MatrixXd& G, const Eigen::VectorXd& f, bool isInequalityConstraint)
    : EqIneqConstraint("Control", isInequalityConstraint)
    , G_(G)
    , f_(f)
{
    if (G.rows() != f.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("G", "f", G, f));
}

void ControlConstraint::reset(const Eigen::MatrixXd& G, const Eigen::VectorXd& f)
{
    if (G.rows() != G_.rows() || G.cols() != G_.cols())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("G", G, G_));
    if (f.rows() != f_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("f", f, f_));

    G_ = G;
    f_ = f;
}

void ControlConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (G_.cols() == ps.uDim) {
        nrConstr_ = static_cast<int>(G_.rows()) * ps.nrUStep;
        A_.resize(nrConstr_, ps.fullUDim);
        b_.resize(nrConstr_);
        A_.setZero();
    } else if (G_.cols() == ps.fullUDim) {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(G_.rows());
        A_ = std::move(G_);
        b_ = std::move(f_);
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSUDim("G", G_, &ps));
    }
}

void ControlConstraint::update(const PreviewSystem& ps)
{
    if (!fullSizeEntry_) {
        auto nrLines = static_cast<int>(G_.rows());
        for (int i = 0; i < ps.nrUStep; ++i) {
            A_.block(i * nrLines, i * ps.uDim, nrLines, ps.uDim) = G_;
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
    if (E.rows() != f.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("E", "f", E, f));
    if (G.rows() != f.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("G", "f", G, f));
}

void MixedConstraint::reset(const Eigen::MatrixXd& E, const Eigen::MatrixXd& G, const Eigen::VectorXd& f)
{
    if (E.rows() != E_.rows() || E.cols() != E_.cols())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("E", E, E_));
    if (E.rows() != E_.rows() || E.cols() != E_.cols())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("G", G, G_));
    if (f.rows() != f_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("f", f, f_));

    E_ = E;
    G_ = G;
    f_ = f;
}

void MixedConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (E_.cols() == ps.xDim) {
        nrConstr_ = static_cast<int>(E_.rows()) * ps.nrUStep;
        A_.resize(nrConstr_, ps.fullUDim);
        b_.resize(nrConstr_);
        A_.setZero();
    } else if (E_.cols() == ps.fullXDim) {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(E_.rows());
        A_.noalias() = E_ * ps.Psi + G_;
        b_.noalias() = f_ - ps.Phi * ps.x0 - ps.xi;
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXUDim("E", "G", E_, G_, &ps));
    }
}

void MixedConstraint::update(const PreviewSystem& ps)
{
    if (!fullSizeEntry_) {
        auto nrLines = static_cast<int>(E_.rows());
        auto uDim = ps.uDim;
        auto xDim = ps.xDim;
        A_.block(0, 0, nrLines, uDim) = G_;
        b_.head(nrLines) = f_ - E_ * ps.x0;
        for (int i = 1; i < ps.nrUStep; ++i) {
            A_.block(i * nrLines, 0, nrLines, uDim) = E_ * ps.Psi.block(i * xDim, 0, xDim, uDim);
            for (int j = 1; j <= i; ++j)
                A_.block(i * nrLines, j * uDim, nrLines, uDim) = A_.block((i - 1) * nrLines, (j - 1) * uDim, nrLines, uDim);

            b_.segment(i * nrLines, nrLines) = f_ - E_ * (ps.Phi.block(i * xDim, 0, xDim, xDim) * ps.x0 + ps.xi.segment(i * xDim, xDim));
        }
    }
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
    if (lower.rows() != upper.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("lower", "upper", lower, upper));

    for (auto line = 0; line < lower_.rows(); ++line) {
        if (lower_(line) != -std::numeric_limits<double>::infinity())
            lowerLines_.push_back(line);
        if (upper_(line) != std::numeric_limits<double>::infinity())
            upperLines_.push_back(line);
    }
}

void TrajectoryBoundConstraint::reset(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
{
    if (lower.rows() != lower_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("lower", lower, lower_));
    if (upper.rows() != upper_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("upper", upper, upper_));

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
    if (lower_.rows() == ps.xDim) {
        nrConstr_ = static_cast<int>((lowerLines_.size() + upperLines_.size())) * ps.nrXStep;
    } else if (lower_.rows() == ps.fullXDim) {
        nrConstr_ = static_cast<int>((lowerLines_.size() + upperLines_.size()));
        fullSizeEntry_ = true;
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSXDim("lower", lower_, &ps));
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
    if (lower.rows() != upper.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("lower", "upper", lower, upper));
}

void ControlBoundConstraint::reset(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper)
{
    if (lower.rows() != lower_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("lower", lower, lower_));
    if (upper.rows() != upper_.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnMat("upper", upper, upper_));

    lower_ = lower;
    upper_ = upper;
}

void ControlBoundConstraint::initializeConstraint(const PreviewSystem& ps)
{
    if (lower_.rows() == ps.uDim) {
        nrConstr_ = ps.fullUDim;
        lb_.resize(nrConstr_);
        ub_.resize(nrConstr_);
    } else if (lower_.rows() == ps.fullUDim) {
        fullSizeEntry_ = true;
        nrConstr_ = static_cast<int>(lower_.rows());
        lb_ = std::move(lower_);
        ub_ = std::move(upper_);
    } else {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnColsOnPSUDim("lower", lower_, &ps));
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