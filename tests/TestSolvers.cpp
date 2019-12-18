/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestSolvers

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

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <numeric>

BOOST_FIXTURE_TEST_CASE(QuadProgTest, Problem)
{
    copra::QuadProgDenseSolver qpQuadProg;

    qpQuadProg.SI_problem(nrvars, nreqs, nrineqs);
    BOOST_REQUIRE(qpQuadProg.SI_solve(Q, c, Aeq, beq, Aineq, bineq, XL, XU));

    BOOST_REQUIRE_EQUAL(qpQuadProg.SI_fail(), 0);
}

#ifdef EIGEN_QLD_FOUND
BOOST_FIXTURE_TEST_CASE(QLDOnQuadProgTest, Problem)
{
    copra::QLDSolver qpQLD;
    copra::QuadProgDenseSolver qpQuadProg;

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

#ifdef EIGEN_LSSOL_FOUND
BOOST_FIXTURE_TEST_CASE(LSSOLOnQuadProgTest, Problem)
{
    copra::QuadProgDenseSolver qpQuadProg;
    copra::LSSOLSolver qpLSSOL;

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

#ifdef EIGEN_GUROBI_FOUND
BOOST_FIXTURE_TEST_CASE(GUROBIOnQuadProgTest, Problem)
{
    copra::QuadProgDenseSolver qpQuadProg;
    copra::GUROBISolver qpGUROBI;

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
