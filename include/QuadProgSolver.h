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

#include "api.h"

#include "SolverInterface.h"
#include <Eigen/Core>
#include <eigen-quadprog/QuadProg.h>

namespace copra {

/**
 * QuadProg solver for dense matrix.
 */
class COPRA_DLLAPI QuadProgDenseSolver : public SolverInterface {
public:
    /**
     * QuadProgDenseSolver default constructor
     */
    QuadProgDenseSolver() = default;

    /**
     * Get information of eventual fail's solver output as define by the
     * solver documentation.
     * \return 0 No problems
     * \return 1 The minimization problem has no solution
     * \return 2 Problems with the decomposition of Q (Is it symmetric positive
     * definite matrix?)
     */
    int SI_fail() const override;

    /**
     * Get the number of needed iteration if available
     * \return The number of iteration
     */
    int SI_iter() const override;

    /**
     * Print an information on the current solver status.
     */
    void SI_inform() const override;

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

    Eigen::QuadProgDense& baseSolver() noexcept { return solver_; }

    QuadProgDenseSolver* clone() const override { return new QuadProgDenseSolver(*this); }

private:
    Eigen::QuadProgDense solver_;
};

} // namespace pc