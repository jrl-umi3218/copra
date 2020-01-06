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
#include <eigen-osqp/OSQP.h>

namespace copra {

/**
 * OSQP solver for dense matrix.
 */

// TODO: Enable sparse matrix
class COPRA_DLLAPI OSQPSolver : public SolverInterface {
public:
    /**
       * QLDSolver default constructor
       */
    OSQPSolver();

    /**
     * Get information of eventual fail's solver output as define by the
     * solver documentation.
     * \return 0 The optimality conditions are satisfied.
     * \return 1 The algorithm has been stopped after too many iterations.
     * \return 2 Termination accuracy insufficient to satisfy convergence
     * criterion.
     * \return 3 Internal inconsistency of QL, division by zero.
     * \return 4 Numerical instability prevents successful termination.
     * \return 5 Length of a working array is too short.
     * \return >100 Constraints are inconsistent and fail=100+ICON, where ICON
     * denotes a constraint causing the conflict.
     */
    int SI_fail() const override;

    void SI_inform() const override;
    int SI_iter() const override;
    int SI_maxIter() const override;
    void SI_maxIter(int maxIter) override;
    void SI_printLevel(int pl) override;
    void SI_feasibilityTolerance(double tol) override;
    bool SI_warmStart() const override;
    void SI_warmStart(bool w) override;

    const Eigen::VectorXd& SI_result() const override;
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;


    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;

    Eigen::OSQP& baseSolver() noexcept { return solver_; }

private:
    Eigen::OSQP solver_;
};

} // namespace pc