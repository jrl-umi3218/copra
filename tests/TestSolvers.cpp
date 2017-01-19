// This file is part of mpc.

// mpc is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with mpc.  If not, see
// <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE TestSolvers

// stl
#include <iostream>
#include <numeric>

// boost
#include <boost/test/unit_test.hpp>

// eigen
#include <Eigen/Core>

// mpc
#include "QuadProgSolver.h"

// optional mpc
#ifdef QLD_SOLVER_FOUND
#include "QLDSolver.h"
#endif
#ifdef LSSOL_SOLVER_FOUND
#include "LSSOLSolver.h"
#endif
#ifdef GUROBI_SOLVER_FOUND
#include "GUROBISolver.h"
#endif

// Test base on scilab qld example:
// https://help.scilab.org/doc/5.5.2/en_US/qld.html
struct Problem {
    Problem()
        : Q(6, 6)
        , Aeq(3, 6)
        , Aineq(2, 6)
        , c(6)
        , beq(3)
        , bineq(2)
        , XL(6)
        , XU(6)
        , nrvars(6)
        , nreqs(3)
        , nrineqs(2)
    {
        Q = Eigen::MatrixXd::Identity(6, 6);
        c << 1, 2, 3, 4, 5, 6;
        Aeq << 1, -1, 1, 0, 3, 1, -1, 0, -3, -4, 5, 6, 2, 5, 3, 0, 1, 0;
        beq << 1, 2, 3;
        Aineq << 0, 1, 0, 1, 2, -1, -1, 0, 2, 1, 1, 0;
        bineq << -1, 2.5;
        XL << -1000, -10000, 0, -1000, -1000, -1000;
        XU << 10000, 100, 1.5, 100, 100, 1000;
    }

    Eigen::MatrixXd Q, Aeq, Aineq;
    Eigen::VectorXd c, beq, bineq, XL, XU;
    int nrvars, nreqs, nrineqs;
};

BOOST_FIXTURE_TEST_CASE(QuadProgTest, Problem)
{
    mpc::QuadProgDenseSolver qpQuadProg;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    BOOST_REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    BOOST_REQUIRE_EQUAL(qpQuadProg.SI_fail(), 0);
}

#ifdef QLD_SOLVER_FOUND
BOOST_FIXTURE_TEST_CASE(QLDOnQuadProgTest, Problem)
{
    mpc::QLDSolver qpQLD;
    mpc::QuadProgDenseSolver qpQuadProg;

    qpQLD.SI_problem(nrvars, nreqs, nrineqs);
    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    BOOST_REQUIRE(qpQLD.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    BOOST_REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQLD = qpQLD.SI_result();
    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    BOOST_CHECK(resQuadProg.isApprox(resQLD));
    BOOST_REQUIRE_EQUAL(qpQLD.SI_fail(), 0);
    BOOST_REQUIRE_EQUAL(qpQuadProg.SI_fail(), 0);
}
#endif

#ifdef LSSOL_SOLVER_FOUND
BOOST_FIXTURE_TEST_CASE(LSSOLOnQuadProgTest, Problem)
{
    mpc::QuadProgDenseSolver qpQuadProg;
    mpc::LSSOLSolver qpLSSOL;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    qpLSSOL.SI_problem(nrvars, nreqs, nrineqs);
    BOOST_REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    BOOST_REQUIRE(qpLSSOL.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resLSSOL = qpLSSOL.SI_result();
    BOOST_CHECK(resLSSOL.isApprox(resQuadProg));
    BOOST_REQUIRE_EQUAL(qpLSSOL.SI_fail(), 0);
}
#endif

#ifdef GUROBI_SOLVER_FOUND
BOOST_FIXTURE_TEST_CASE(GUROBIOnQuadProgTest, Problem)
{
    mpc::QuadProgDenseSolver qpQuadProg;
    mpc::GUROBISolver qpGUROBI;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    qpGUROBI.SI_problem(nrvars, nreqs, nrineqs);
    BOOST_REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    BOOST_REQUIRE(qpGUROBI.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resGUROBI = qpGUROBI.SI_result();
    BOOST_CHECK(resGUROBI.isApprox(resQuadProg, 1e-6));
    BOOST_REQUIRE_EQUAL(qpGUROBI.SI_fail(), GRB_OPTIMAL);
}
#endif