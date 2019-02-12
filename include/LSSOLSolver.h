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
#include <eigen-lssol/LSSOL_QP.h>

namespace copra {

/**
 * LSSOLSolver solver for both dense matrix.
 */
class COPRA_DLLAPI LSSOLSolver : public SolverInterface {
public:
    /**
     * LSSOLSolver default constructor
     */
    LSSOLSolver() = default;

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
    /**
     * Select the print level of the solver if available
     * \param pl The print level.
     * \param pl =0  : Nothing is printed
     * \param pl =1  : The final solution is printed
     * \param pl =5  : One line of output for each iteration
     * \param pl =10 : One line of output for each iteration and the final solution is printed
     * \param pl >20 : At each iteration, the Lagrangian multipliers, the variables x,
     * the constraints values Cx and the constraints status.
     * \param pl >30 : For an understanding of this, see the official documentation.
     */
    void SI_printLevel(int pl) override;

    void SI_inform() const override;
    int SI_iter() const override;

    int SI_maxIter() const override;
    void SI_maxIter(int maxIter) override;

    double SI_feasibilityTolerance() const override;
    void SI_feasibilityTolerance(double tol) override;

    bool SI_warmStart() const override;
    void SI_warmStart(bool w) override;

    const Eigen::VectorXd& SI_result() const override;
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;
    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& c,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;

    Eigen::LSSOL_QP& baseSolver() noexcept { return solver_; }

    LSSOLSolver* clone() const override { return new LSSOLSolver(*this); }

private:
    Eigen::LSSOL_QP solver_;
    Eigen::MatrixXd Q_;
    Eigen::MatrixXd A_;
    Eigen::VectorXd bl_;
    Eigen::VectorXd bu_;
};

} // namespace pc
