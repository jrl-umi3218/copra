#include "SolverInterface.h"

// std
#include <utility>

namespace pc
{

/*
 * SolverInterface
 */
int SolverInterface::SI_fail() const
{
    return 0;
}

const Eigen::VectorXd &SolverInterface::SI_result() const
{
    return std::move(Eigen::VectorXd());
}

void SolverInterface::SI_problem(int /* nrVar */, int /* nrEq */, int /* nrInEq */)
{
}

bool SolverInterface::SI_solve(const Eigen::MatrixXd & /* Q */, const Eigen::VectorXd & /* C */,
                               const Eigen::MatrixXd & /* Aeq */, const Eigen::VectorXd & /* Beq */,
                               const Eigen::MatrixXd & /* Aineq */, const Eigen::VectorXd & /* Bineq */,
                               const Eigen::VectorXd & /* XL */, const Eigen::VectorXd & /* XU */)
{
    return false;
}

} // namespace pc