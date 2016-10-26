//This file is part of ModelPreviewController.

//ModelPreviewController is free software: you can redistribute it and/or modify
//it under the terms of the GNU Lesser General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//ModelPreviewController is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU Lesser General Public License for more details.
//
//You should have received a copy of the GNU Lesser General Public License
//along with ModelPreviewController.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <Eigen/Core>
#include <string>

namespace mpc
{

//forward declaration
struct PreviewSystem;

/**
 * Constraint to add to the system.
 * This is a pure virtual class
 */
class Constraint
{
  public:
    /**
     * Constructor of a constraint.
     * @param E The inequality constraint matrix.
     * @param f The inequality constraint vector.
     */
    Constraint(const Eigen::MatrixXd &E, const Eigen::VectorXd &f);

    /**
     * Declare virtual desturctor
     */
    virtual ~Constraint() = default;

    /**
     * Initialization of the constraint.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void initializeConstraint(const PreviewSystem &ps) = 0;
    /**
     * Update the constraint.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void update(const PreviewSystem &ps) = 0;

    /**
     * Function that return the number of constrains.
     * @return The number of constrains.
     */
    virtual std::string name() const noexcept = 0;

    /**
     * Function that return the number of constrains.
     * @return The number of constrains.
     */
    int nrConstr() noexcept
    {
        return nrConstr_;
    }
    /**
     * Function that return the QP form-like inequality matrix.
     * @return The inequality matrix.
     */
    const Eigen::MatrixXd &A() noexcept
    {
        return A_;
    }
    /**
     * Function that return the QP form-like inequality vector.
     * @return The inequality vector.
     */
    const Eigen::VectorXd &b() noexcept
    {
        return b_;
    }

  protected:
    std::string name_;
    int nrConstr_;
    Eigen::MatrixXd E_, A_;
    Eigen::VectorXd f_, b_;
};

class TrajectoryConstraint final : public Constraint
{
  public:
    /**
     * Constructor of a constraint.
     * Constraint of type \f$EX <= f\f$ 
     * @param E The inequality constraint matrix.
     * @param f The inequality constraint vector.
     */
    using Constraint::Constraint; //Inherits base class constructor
    std::string name() const noexcept override;
    void initializeConstraint(const PreviewSystem &ps) override;
    void update(const PreviewSystem &ps) override;
};

class ControlConstraint final : public Constraint
{
  public:
    /**
     * Constructor of a constraint.
     * Constraint of type \f$EU <= f\f$ 
     * @param E The inequality constraint matrix.
     * @param f The inequality constraint vector.
     */
    using Constraint::Constraint; //Inherits base class constructor
    std::string name() const noexcept override;
    void initializeConstraint(const PreviewSystem &ps) override;
    void update(const PreviewSystem &ps) override;
};

} // namespace mpc