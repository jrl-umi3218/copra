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

#pragma once

#include <Eigen/Core>
#include <string>

namespace mpc {

//forward declaration
struct PreviewSystem;

enum class ConstraintFlag {
    EqualityConstraint,
    InequalityConstraint,
    BoundConstraint
};

/**
 * Constraint to add to the system.
 * This is a pure virtual class
 */
class Constraint {
public:
    /**
     * Constructor of a constraint.
     */
    Constraint(const std::string& name);

    /**
     * Declare virtual desturctor
     */
    virtual ~Constraint() = default;

    /**
     * Initialization of the constraint.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void initializeConstraint(const PreviewSystem& ps) = 0;

    /**
     * Update the constraint.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void update(const PreviewSystem& ps) = 0;

    /**
     * Get the type of the constraint
     * @return The type of the constraint @see ConstraintFlag
     */
    virtual ConstraintFlag constraintType() = 0;

    /**
     * Function that return the number of constrains.
     * @return The number of constrains.
     */
    const std::string& name() const noexcept
    {
        return name_;
    }

    /**
     * Function that return the number of constrains.
     * @return The number of constrains.
     */
    int nrConstr() noexcept
    {
        return nrConstr_;
    }

protected:
    std::string name_;
    int nrConstr_;
};

class EqIneqConstraint : public Constraint {
public:
  EqIneqConstraint(const std::string &name, bool isInequalityConstraint);

  const Eigen::MatrixXd &A()
  {
      return A_;
    }

    const Eigen::VectorXd& b()
    {
        return b_;
    }

protected:
    Eigen::MatrixXd A_;
    Eigen::VectorXd b_;
    bool isIneq_;
};

class TrajectoryConstraint final : public EqIneqConstraint {
public:
    TrajectoryConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint = true);
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    ConstraintFlag constraintType() override;

private:
    Eigen::MatrixXd E_;
    Eigen::VectorXd f_;
};

class ControlConstraint final : public EqIneqConstraint {
public:
    ControlConstraint(const Eigen::MatrixXd& E, const Eigen::VectorXd& f, bool isInequalityConstraint = true);
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    ConstraintFlag constraintType() override;

private:
    Eigen::MatrixXd E_;
    Eigen::VectorXd f_;
    bool isIneq_;
};

class TrajectoryBoundConstraint final : public EqIneqConstraint {
public:
    TrajectoryBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper);
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    ConstraintFlag constraintType() override;

private:
    Eigen::VectorXd lower_, upper_;
};

class ControlBoundConstraint final : public Constraint {
public:
    ControlBoundConstraint(const Eigen::VectorXd& lower, const Eigen::VectorXd& upper);
    void initializeConstraint(const PreviewSystem& ps) override;
    void update(const PreviewSystem& ps) override;
    ConstraintFlag constraintType() override;

    const Eigen::VectorXd& lower()
    {
        return lb_;
    }

    const Eigen::VectorXd& upper()
    {
        return ub_;
    }

private:
    Eigen::VectorXd lower_, upper_, lb_, ub_;
};

} // namespace mpc