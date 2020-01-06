/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include <Eigen/Core>

namespace copra {

/**
 * An interface to the quadratic solvers.
 * This interface is more like a pseudo-interface (it is not an abstract class).
 * This class allows to have a base class for all solvers.
 * It provides all the necessary functions for using a qp solver.
 * In case of a QP does not have a corresponding function it sends a warning.
 */
class COPRA_DLLAPI SolverInterface {
public:
    /**
     * Default virtual destructor.
     */
    virtual ~SolverInterface() = default;
    /**
     * Get information of eventual fail's solver output as define by the
     * solver documentation.
     * \return The fail number.
     */
    virtual int SI_fail() const = 0;

    /**
     * Print an information on the current solver status.
     */
    virtual void SI_inform() const = 0;

    /**
     * Get the number of needed iteration if available.
     * \return The number of iteration
     */
    virtual int SI_iter() const;

    /**
     * Get max number of iteration if available.
     * \return Max number of iteration
     */
    virtual int SI_maxIter() const;

    /**
     * Set max number of iteration if available.
     */
    virtual void SI_maxIter(int maxIter);

    /**
     * Select the print level of the solver if available.
     * \param pl The print level
     */
    virtual void SI_printLevel(int pl);

    /**
     * Get the maximum error tolerance of the solution if available.
     */
    virtual double SI_feasibilityTolerance() const;

    /**
     * Set the maximum error tolerance of the solution if available.
     * \param tol The error tolerance
     */
    virtual void SI_feasibilityTolerance(double tol);

    /**
     * Get the warm start status of the solver if available.
     * \return True for warm start, False for cold start
     */
    virtual bool SI_warmStart() const;

    /**
     * Set the warm start status of the solver if available.
     * \param w True for warm start, False for cold start
     */
    virtual void SI_warmStart(bool w);

    /**
     * Get the solver's solution.
     * \return The qp solver result.
     */
    virtual const Eigen::VectorXd& SI_result() const = 0;

    /**
     * Initialize the variables of the problem to solve.
     * \param nrVar The number of decision variables
     * \param nrEq The number of equality constrains
     * \param nrInEq The number of inequality constrains
     * \return The qp solver result.
     */
    virtual void SI_problem(int nrVar, int nrEq, int nrInEq) = 0;

    /**
     * Solve the problem.
     * Solve the system:
     \f[
      \left\{
        \begin{array}{cl}
          \min\limits_x & \frac{1}{2} x^T Q x + c^T x \\
                        & A_{eq} x = b_{eq} \\
                        & A_{ineq} x \leq b_{ineq} \\
                        & XL \leq x \leq XU
        \end{array}
      \right.
     \f]
     * \param Q An N-by-N symmetric positive definite dense matrix.
     * \param c An N-by-1 dense vector.
     * \param Aeq Equality constrain. An M-by-N dense matrix.
     * \param beq Equality constrain. An M-by-1 dense vector.
     * \param Aineq Inequality constrain. An P-by-N dense matrix.
     * \param bineq Inequality constrain. An P-by-1 dense vector.
     * \param XL Lower bound of the decision vector. An N-by-1 dense vector.
     * \param XU Upper bound of the decision vector. An N-by-1 dense vector.
     * \return The qp solver result.
     */
    virtual bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU)
        = 0;
};

} // namespace pc
