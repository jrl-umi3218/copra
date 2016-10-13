#pragma once

#include <memory>
#include <Eigen/Core>
#include <eigen-qld/QLD.h>
#include "SolverInterface.h"

namespace pc
{

/**
 * QLD solver for both dense matrix.
 */
class QLDSolver : public SolverInterface //TODO: Enable sparse matrix
{
  public:
    /**
	 * QLDSolver default constructor
	 */
    QLDSolver();

    /**
	 * Get information of eventual fail's solver output as define by the solver documentation.
	 * @return 0 The optimality conditions are satisfied. 
	 * @return 1 The algorithm has been stopped after too many iterations. 
	 * @return 2 Termination accuracy insufficient to satisfy convergence criterion. 
	 * @return 3 Internal inconsistency of QL, division by zero. 
	 * @return 4 Numerical instability prevents successful termination.
	 * @return 5 Length of a working array is too short. 
	 * @return >100 Constraints are inconsistent and fail=100+ICON, where ICON denotes a constraint causing the conflict. 
	 */
    int SI_fail() const override;

    /**
	 * Get the solver's solution.
	 * @return The qp solver result.
	 */
    const Eigen::VectorXd &SI_result() const override;

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
    bool SI_solve(const Eigen::MatrixXd &Q, const Eigen::VectorXd &C,
                  const Eigen::MatrixXd &Aeq, const Eigen::VectorXd &Beq,
                  const Eigen::MatrixXd &Aineq, const Eigen::VectorXd &Bineq,
                  const Eigen::VectorXd &XL, const Eigen::VectorXd &XU) override;

  private:
    std::unique_ptr<Eigen::QLD> solver_;
};

} // namespace pc