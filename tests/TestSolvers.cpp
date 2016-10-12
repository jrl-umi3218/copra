#define BOOST_TEST_MODULE Hello 
#include <boost/test/unit_test.hpp>
#include "solvers.h"
#include <Eigen/Core>
#include <numeric>

BOOST_AUTO_TEST_CASE(SolverTest)
{
    pc::QLDSolver qp;

    Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(3, 3);
    Eigen::VectorXd c = Eigen::VectorXd::Zero(3);
    Eigen::MatrixXd Aineq = -Eigen::MatrixXd::Ones(0, 3);
    Eigen::VectorXd bineq = Eigen::VectorXd::Zero(0);
    Eigen::MatrixXd Aeq = Eigen::MatrixXd::Zero(0, 3);
    Eigen::VectorXd beq = Eigen::VectorXd::Zero(0);
    Eigen::VectorXd XL = Eigen::VectorXd::Constant(3, 1);
    Eigen::VectorXd XU = Eigen::VectorXd::Constant(3, std::numeric_limits<double>::infinity());

    qp.SI_problem(3, 0, 3);
    BOOST_REQUIRE(qp.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd res = qp.SI_result();
    BOOST_CHECK(res.isApproxToConstant(1));
}