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

#include "QuadProgSolver.h"
#include "systems.h"
#ifdef EIGEN_QLD_FOUND
#include "QLDSolver.h"
#endif
#ifdef EIGEN_LSSOL_FOUND
#include "LSSOLSolver.h"
#endif
#ifdef EIGEN_GUROBI_FOUND
#include "GUROBISolver.h"
#endif
#ifdef EIGEN_OSQP_FOUND
#include "OSQPSolver.h"
#endif

#include <catch2/catch.hpp>
#include <Eigen/Core>
#include <iostream>
#include <numeric>

#define REQUIRE_EQUAL(a, b) REQUIRE((a) == (b))

TEST_CASE_METHOD(Problem, "QuadProgTest", "[dense]")
{
    copra::QuadProgDenseSolver qpQuadProg;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    REQUIRE_EQUAL(qpQuadProg.SI_fail(), 0);
}

#ifdef EIGEN_QLD_FOUND
TEST_CASE_METHOD(Problem, "QLDOnQuadProgTest", "[dense]")
{
    copra::QLDSolver qpQLD;
    copra::QuadProgDenseSolver qpQuadProg;

    qpQLD.SI_problem(nrvars, nreqs, nrineqs);
    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    REQUIRE(qpQLD.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQLD = qpQLD.SI_result();
    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    CHECK(resQuadProg.isApprox(resQLD));
    REQUIRE_EQUAL(qpQLD.SI_fail(), 0);
    REQUIRE_EQUAL(qpQuadProg.SI_fail(), 0);
}
#endif

#ifdef EIGEN_LSSOL_FOUND
TEST_CASE_METHOD(Problem, "LSSOLOnQuadProgTest", "[dense]")
{
    copra::QuadProgDenseSolver qpQuadProg;
    copra::LSSOLSolver qpLSSOL;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    qpLSSOL.SI_problem(nrvars, nreqs, nrineqs);
    REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    REQUIRE(qpLSSOL.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resLSSOL = qpLSSOL.SI_result();
    CHECK(resLSSOL.isApprox(resQuadProg));
    REQUIRE_EQUAL(qpLSSOL.SI_fail(), 0);
}
#endif

#ifdef EIGEN_GUROBI_FOUND
TEST_CASE_METHOD(Problem, "GUROBIOnQuadProgTest", "[dense]")
{
    copra::QuadProgDenseSolver qpQuadProg;
    copra::GUROBISolver qpGUROBI;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    qpGUROBI.SI_problem(nrvars, nreqs, nrineqs);
    REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    REQUIRE(qpGUROBI.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resGUROBI = qpGUROBI.SI_result();
    CHECK(resGUROBI.isApprox(resQuadProg, 1e-6));
    REQUIRE_EQUAL(qpGUROBI.SI_fail(), GRB_OPTIMAL);
}
#endif

#ifdef EIGEN_OSQP_FOUND
TEST_CASE_METHOD(Problem, "OSQPOnQuadProgTest", "[dense]")
{
    copra::QuadProgDenseSolver qpQuadProg;
    copra::OSQPSolver osqp;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    osqp.SI_problem(nrvars, nreqs, nrineqs);
    REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));
    REQUIRE(osqp.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    Eigen::VectorXd resQuadProg = qpQuadProg.SI_result();
    Eigen::VectorXd resOSQP = osqp.SI_result();
    CHECK(resOSQP.isApprox(resQuadProg, 1e-6));
    REQUIRE_EQUAL(osqp.SI_fail(), 0);
}
#endif
