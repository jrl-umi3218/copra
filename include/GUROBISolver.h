// This file is part of copra.

// copra is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// copra is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with copra.  If not, see
// <http://www.gnu.org/licenses/>.

#pragma once

// stl
#include <memory>

// eigen
#include <Eigen/Core>
#include <eigen-gurobi/Gurobi.h>

// copra
#include "SolverInterface.h"
#include "config.hh"

namespace copra {

/**
 * GUROBISolver solver for both dense matrix.
 */
class COPRA_DLLAPI GUROBISolver : public SolverInterface // TODO: Enable sparse matrix
{
public:
    /**
     * GUROBISolver default constructor
     */
    GUROBISolver();

    /**
     * Get information of eventual fail's solver output as define by the
     * solver documentation.
     * See the solver documentation (Need to add the doc here)
     */
    int SI_fail() const override;

    /**
	 * Print details on the status of the solver as described in the docs.
	 */
    void SI_inform() const override;

    /**
	 * Get the number of needed iteration if available
	 * \return The number of iteration
	 */
    int SI_iter() const override;

    /**
     * Display or hide (default) the output of the solver.
     * There is no printing level. It shows everything or nothing.
     * \param pl 0 (hide) or else (display)
     */
    void SI_printLevel(int pl) override;

    /**
     * Set the maximum error tolerance of the solution
     * \param tol The error tolerance
     */
    void SI_feasibilityTolerance(double tol) override;

    /**
     * Get the warm start status of the solver
     * \return True for warm start, False for cold start
     */
    bool SI_warmStart() const override;

    /**
     * Set the warm start status of the solver
     * \param w True for warm start, False for cold start
     */
    void SI_warmStart(bool w) override;

    /**
     * Get the solver's solution.
     * \return The qp solver result.
     */
    const Eigen::VectorXd& SI_result() const override;

    /**
     * Initialize the variables of the problem to solve.
     * \see SolverInterface::SI_problem()
     * \return The qp solver result.
     */
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;

    /**
     * Solve the problem.
     * \see SolverInterface::SI_solve()
     * \return The qp solver result.
     */
    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;

private:
    std::unique_ptr<Eigen::GurobiDense> solver_;
};

} // namespace pc