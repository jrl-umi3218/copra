#pragma once

#include <Eigen/Core>

namespace mpc
{

//Forward declare
class PreviewSystem;

/**
 * Constrain to add to the system.
 * This is a pure virtual class
 */
class Constrain
{
  public:
    /**
     * Constructor of a constrain.
     * @param E The inequality constrain matrix.
     * @param f The inequality constrain vector.
     */
    Constrain(const Eigen::MatrixXd &E, const Eigen::VectorXd &f);

    /**
     * Declare virtual desturctor
     */
    virtual ~Constrain() = default;

    /**
     * Initialization of the constrain.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void initializeConstrain(const PreviewSystem &ps) = 0;
    /**
     * Update the constrain.
     * @param ps An unique pointer to a PreviewSystem.
     */
    virtual void update(const PreviewSystem &ps) = 0;

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
    int nrConstr_;
    Eigen::MatrixXd E_, A_;
    Eigen::VectorXd f_, b_;
};

class TrajectoryConstrain final : public Constrain
{
  public:
    /**
     * Constructor of a constrain.
     * Constrain of type \f$EX <= f\f$ 
     * @param E The inequality constrain matrix.
     * @param f The inequality constrain vector.
     */
    using Constrain::Constrain; //Inherits base class constructor
    void initializeConstrain(const PreviewSystem &ps) override;
    void update(const PreviewSystem &ps) override;
};

class ControlConstrain final : public Constrain
{
  public:
    /**
     * Constructor of a constrain.
     * Constrain of type \f$EU <= f\f$ 
     * @param E The inequality constrain matrix.
     * @param f The inequality constrain vector.
     */
    using Constrain::Constrain; //Inherits base class constructor
    void initializeConstrain(const PreviewSystem &ps) override;
    void update(const PreviewSystem &ps) override;
};

} // namespace mpc