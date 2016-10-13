#pragma once

#include "SolverInterface.h"
#include <Eigen/Core>
#include <eigen-lssol/LSSOL.h>
#include <memory>

namespace pc {

/**
 * LSSOLSolver solver for both dense matrix.
 */
class LSSOLSolver : public SolverInterface // TODO: Enable sparse matrix
{
public:
    /**
       * LSSOLSolver default constructor
       */
    LSSOLSolver();

    /**
       * Get information of eventual fail's solver output as define by the
   * solver documentation.
       * @return 0 The optimality conditions are satisfied.
       * @return 1 The algorithm has been stopped after too many iterations.
       * @return 2 Termination accuracy insufficient to satisfy convergence
   * criterion.
       * @return 3 Internal inconsistency of QL, division by zero.
       * @return 4 Numerical instability prevents successful termination.
       * @return 5 Length of a working array is too short.
       * @return >100 Constraints are inconsistent and fail=100+ICON, where ICON
   * denotes a constraint causing the conflict.
       */
    int SI_fail() const override;
    int SI_iter() const override;
    void SI_inform() const override;
    void SI_printLevel(int pl) const override;
    void SI_tol(double tol) const override;
    bool SI_warmStart() const override;
    void SI_warmStart(bool w) override;

    /**
       * Get the solver's solution.
       * @return The qp solver result.
       */
    const Eigen::VectorXd& SI_result() const override;

    /**
       * Initialize the variables of the problem to solve.
       * @see SolverInterface::SI_problem()
       * @return The qp solver result.
       */
    void SI_problem(int nrVar, int nrEq, int nrInEq) override;

    /**
       * Solve the problem.
       * @see SolverInterface::SI_solve()
       * @return The qp solver result.
       */
    bool SI_solve(const Eigen::MatrixXd& Q, const Eigen::VectorXd& C,
        const Eigen::MatrixXd& Aeq, const Eigen::VectorXd& Beq,
        const Eigen::MatrixXd& Aineq, const Eigen::VectorXd& Bineq,
        const Eigen::VectorXd& XL, const Eigen::VectorXd& XU) override;

private:
    std::unique_ptr<Eigen::StdLSSOL> solver_;
};

} // namespace pc