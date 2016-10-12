#include "solver.h"
#include <iostream>
#include <Eigen/Core>

int main()
{
    QuadProgDenseSolver qp();

    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(3, 3);
    Eigen::VectorXd c = Eigen::VectorXd::Zero(3);
    Eigen::MatrixXd Aineq = -Eigen::MatrixXd::Ones(1, 3);
    Eigen::VectorXd bineq = 1000*Eigen::VectorXd::Ones(3);
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(0, 3);
    Eigen::VectorXd beq = Eigen::VectorXd::Zero(0);

    qp.SI_problem(3, 0, 3);
    bool success = qp.SI_solve(Q, c, Aeq, beq, Aineq, bineq);
    Eigen::VectorXd res = qp.si_result();

    for(auto i = 0; i < res.rows(); ++i)
        std::cout << res(i) << std::endl;

    return 0;
}