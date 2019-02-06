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

#define BOOST_TEST_MODULE TestLMPC

#include "systems.h"
#include "tools.h"
#include <Eigen/Core>
#include <algorithm>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include "LMPC.h"
#include "PreviewSystem.h"
#include "constraints.h"
#include "costFunctions.h"
#include "QuadProgSolver.h"
#ifdef EIGEN_QLD_FOUND
#include "QLDSolver.h"
#endif
#ifdef EIGEN_LSSOL_FOUND
#include "LSSOLSolver.h"
#endif
#ifdef EIGEN_GUROBI_FOUND
#include "GUROBISolver.h"
#endif
#include <memory>
#include <numeric>
#include <vector>

/********************************************************************************************************
 *                               Check Bound constraint                                                 *
 ********************************************************************************************************/

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_BOUND_CONSTRAINTS, BoundedSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TargetCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
        BOOST_REQUIRE_LE(velTraj.maxCoeff(), xUpper(1) + 1e-6);
        BOOST_REQUIRE_LE(control.maxCoeff(), uUpper(0) + 1e-6);
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_BOUND_CONSTRAINTS, BoundedSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
        BOOST_REQUIRE_LE(velTraj.maxCoeff(), xUpper(1) + 1e-6);
        BOOST_REQUIRE_LE(control.maxCoeff(), uUpper(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_MIXED_COST_WITH_BOUND_CONSTRAINTS, BoundedSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd); // min(||X - Xt||^2)
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud); // min(||U - Ut||^2)
    auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    BOOST_REQUIRE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check result
    BOOST_CHECK_SMALL(xd(1) - velTraj.tail(3)(0), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
    BOOST_REQUIRE_LE(velTraj.maxCoeff(), xUpper(1) + 1e-6);
    BOOST_REQUIRE_LE(control.maxCoeff(), uUpper(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
}

/********************************************************************************************************
 *                            Check inequality constraint                                               *
 ********************************************************************************************************/

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_INEQUALITY_CONSTRAINTS, IneqSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TargetCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f);
    auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
        BOOST_REQUIRE_LE(velTraj.maxCoeff(), f(0) + 1e-6);
        BOOST_REQUIRE_LE(control.maxCoeff(), h(0) + 1e-6);
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_INEQUALITY_CONSTRAINTS, IneqSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f);
    auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
        BOOST_REQUIRE_LE(velTraj.maxCoeff(), f(0) + 1e-6);
        BOOST_REQUIRE_LE(control.maxCoeff(), h(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_MIXED_COST_WITH_INEQUALITY_CONSTRAINTS, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd);
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f);
    auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    BOOST_REQUIRE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check result
    BOOST_CHECK_SMALL(xd(1) - velTraj.tail(3)(0), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
    BOOST_REQUIRE_LE(velTraj.maxCoeff(), f(0) + 1e-6);
    BOOST_REQUIRE_LE(control.maxCoeff(), h(0) + 1e-6); // QuadProg allows to exceeds the constrain of a small amount.
}

/********************************************************************************************************
 *                               Check Mixed constraint                                                 *
 ********************************************************************************************************/

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_MIXED_CONSTRAINTS, MixedSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TargetCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, f);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
        for (int i = 0; i < nbStep; ++i) {
            auto res = E * fullTraj.segment(i * E.cols(), E.cols()) + G * control.segment(i * G.cols(), G.cols());
            if (!(res(0) <= f(0) + 1e-6))
                BOOST_ERROR("Mixed constraint violated!");
        }
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_MIXED_CONSTRAINTS, MixedSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, f);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
        for (int i = 0; i < nbStep; ++i) {
            auto res = E * fullTraj.segment(i * E.cols(), E.cols()) + G * control.segment(i * G.cols(), G.cols());
            if (!(res(0) <= f(0) + 1e-6))
                BOOST_ERROR("Mixed constraint violated!");
        }
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_MIXED_COST_WITH_MIXED_CONSTRAINTS, MixedSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd);
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud);
    auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, f);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    BOOST_REQUIRE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check result
    BOOST_CHECK_SMALL(xd(1) - velTraj.tail(3)(0), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
    for (int i = 0; i < nbStep; ++i) {
        auto res = E * fullTraj.segment(i * E.cols(), E.cols()) + G * control.segment(i * G.cols(), G.cols());
        if (!(res(0) <= f(0) + 1e-6))
            BOOST_ERROR("Mixed constraint violated!");
    }
}

/********************************************************************************************************
 *                              Check Equality constraint                                               *
 ********************************************************************************************************/

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_EQUALITY_CONSTRAINTS, EqSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TargetCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0) + 1e-6);
        BOOST_REQUIRE_LE(velTraj.maxCoeff(), f(0) + 1e-6);
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    solver->SI_feasibilityTolerance(1e-6);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_EQUALITY_CONSTRAINTS, EqSystem)
{
    tools::SolverTimers sTimers;

    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::TrajectoryCost>(M, xd);
    auto uCost = std::make_shared<copra::ControlCost>(N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    auto pcCheck = [&](const std::string& solverName, copra::SolverFlag sFlag, std::unique_ptr<copra::SolverInterface>&& solver = nullptr) {
        BOOST_TEST_MESSAGE(solverName);
        if (solver)
            controller.useSolver(std::move(solver));
        else
            controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        sTimers.st.emplace_back(solverName, controller.solveTime() * 1e3);
        sTimers.ct.emplace_back(solverName, controller.solveAndBuildTime() * 1e3);
        sTimers.bt.emplace_back(solverName, (controller.solveAndBuildTime() - controller.solveTime()) * 1e3);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i) {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0) + 1e-6);
        BOOST_REQUIRE_LE(velTraj.maxCoeff(), f(0) + 1e-6);
    };

    pcCheck("Default (QuadProgDense)", copra::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    auto solver = copra::solverFactory(copra::SolverFlag::LSSOL);
    solver->SI_maxIter(130);
    solver->SI_feasibilityTolerance(1e-6);
    pcCheck("LSSOL", copra::SolverFlag::LSSOL, std::move(solver));
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", copra::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", copra::SolverFlag::GUROBIDense);
#endif

    tools::printSortedTimers(sTimers);
}

BOOST_FIXTURE_TEST_CASE(MPC_MIXED_COST_WITH_EQUALITY_CONSTRAINTS, EqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto xCost = std::make_shared<copra::MixedCost>(M, Eigen::MatrixXd::Zero(2, 1), xd);
    auto uCost = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Zero(1, 2), N, ud);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    BOOST_REQUIRE(controller.solve());

    Eigen::VectorXd fullTraj = controller.trajectory();
    auto trajLen = fullTraj.rows() / 2;
    Eigen::VectorXd posTraj(trajLen);
    Eigen::VectorXd velTraj(trajLen);
    for (auto i = 0; i < trajLen; ++i) {
        posTraj(i) = fullTraj(2 * i);
        velTraj(i) = fullTraj(2 * i + 1);
    }
    Eigen::VectorXd control = controller.control();

    // Check result
    BOOST_CHECK_SMALL(xd(1) - velTraj.tail(3)(0), 0.001); // Check X_{N-1} for mixed cost because X_N is not evaluated.

    // Check constrains
    BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0) + 1e-6);
    BOOST_REQUIRE_LE(velTraj.maxCoeff(), f(0) + 1e-6);
}

/********************************************************************************************************
 *                                   Check Autospan                                                     *
 ********************************************************************************************************/

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_BOUND_CONSTRAINT, BoundedSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::VectorXd& xLower, const Eigen::VectorXd& xUpper, const Eigen::VectorXd& uLower, const Eigen::VectorXd& uUpper) {
        auto trajConstr = std::make_shared<copra::TrajectoryBoundConstraint>(xLower, xUpper);
        trajConstr->autoSpan();

        auto contConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
        contConstr->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addConstraint(trajConstr));
        BOOST_REQUIRE_NO_THROW(controller.addConstraint(contConstr));
    };

    auto fullxLower = tools::spanVector(xLower, nbXStep);
    auto fullxUpper = tools::spanVector(xUpper, nbXStep);
    auto fulluLower = tools::spanVector(uLower, nbStep);
    auto fulluUpper = tools::spanVector(uUpper, nbStep);

    checkSpan(xLower, xUpper, uLower, uUpper);
    checkSpan(fullxLower, xUpper, fulluLower, uUpper);
    checkSpan(xLower, fullxUpper, uLower, fulluUpper);
    checkSpan(fullxLower, fullxUpper, fulluLower, fulluUpper);
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_INEQUALITY_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::MatrixXd& E, const Eigen::VectorXd& f, const Eigen::MatrixXd& G, const Eigen::VectorXd& h) {
        auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f);
        trajConstr->autoSpan();

        auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);
        contConstr->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addConstraint(trajConstr));
        BOOST_REQUIRE_NO_THROW(controller.addConstraint(contConstr));
    };

    auto fullE = tools::spanMatrix(E, nbXStep);
    auto fullf = tools::spanVector(f, nbXStep);
    auto fullG = tools::spanMatrix(G, nbStep);
    auto fullh = tools::spanVector(h, nbStep);

    checkSpan(E, f, G, h);
    checkSpan(fullE, f, fullG, h);
    checkSpan(E, fullf, G, fullh);
    checkSpan(fullE, fullf, fullG, fullh);
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_MIXED_CONSTRAINT, MixedSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto checkSpan = [&](const Eigen::MatrixXd& E, const Eigen::MatrixXd& G, const Eigen::VectorXd& f) {
        auto mixedConstr = std::make_shared<copra::MixedConstraint>(E, G, f);
        mixedConstr->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addConstraint(mixedConstr));
    };

    auto fullE = tools::spanMatrix(E, nbStep, 1);
    auto fullG = tools::spanMatrix(G, nbStep);
    auto fullf = tools::spanVector(f, nbStep);

    checkSpan(E, G, f);
    checkSpan(fullE, G, f);
    checkSpan(fullE, fullG, f);
    checkSpan(fullE, G, fullf);
    checkSpan(E, fullG, f);
    checkSpan(fullE, fullG, f);
    checkSpan(E, fullG, fullf);
    checkSpan(E, G, fullf);
    checkSpan(fullE, G, fullf);
    checkSpan(E, fullG, fullf);
    checkSpan(fullE, fullG, fullf);
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_TRAJECTORY_COST, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::MatrixXd& M, const Eigen::VectorXd& p, const Eigen::VectorXd& weights) {
        auto cost = std::make_shared<copra::TrajectoryCost>(M, p);
        cost->weights(weights);
        cost->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addCost(cost));
    };

    auto fullM = tools::spanMatrix(M, nbXStep);
    auto fullxd = tools::spanVector(xd, nbXStep);

    checkSpan(M, xd, wx);
    checkSpan(M, fullxd, wx);
    checkSpan(fullM, xd, wx);
    checkSpan(fullM, fullxd, wx);
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_CONTROL_COST, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto checkSpan = [&](const Eigen::MatrixXd& M, const Eigen::VectorXd& p, const Eigen::VectorXd& weights) {
        auto cost = std::make_shared<copra::ControlCost>(M, p);
        cost->weights(weights);
        cost->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addCost(cost));
    };

    auto fullN = tools::spanMatrix(N, nbStep);
    auto fullud = tools::spanVector(ud, nbStep);

    checkSpan(N, ud, wu);
    checkSpan(N, fullud, wu);
    checkSpan(fullN, ud, wu);
    checkSpan(fullN, fullud, wu);
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_MIXED_COST, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto checkSpan = [&](const Eigen::MatrixXd& M, const Eigen::MatrixXd& N, const Eigen::VectorXd& p, const Eigen::VectorXd& weights) {
        auto cost = std::make_shared<copra::MixedCost>(M, N, p);
        cost->weights(weights);
        cost->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addCost(cost));
    };

    auto MVec = std::vector<Eigen::MatrixXd>();
    MVec.push_back(M);
    MVec.push_back(tools::spanMatrix(M, nbStep, 1));
    auto nnVec = std::vector<Eigen::MatrixXd>();
    nnVec.push_back(Eigen::MatrixXd::Ones(2, 1));
    nnVec.push_back(tools::spanMatrix(Eigen::MatrixXd::Ones(2, 1), nbStep));
    auto xdVec = std::vector<Eigen::VectorXd>();
    xdVec.push_back(xd);
    xdVec.push_back(tools::spanVector(xd, nbStep));

    for (auto& i : MVec)
        for (auto& j : nnVec)
            for (auto& k : xdVec)
                checkSpan(i, j, k, wx);
}

/********************************************************************************************************
 *                                Check Error Messages                                                  *
 ********************************************************************************************************/

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_PREVIEW_SYSTEM, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    BOOST_REQUIRE_THROW(ps->system(Eigen::MatrixXd::Ones(5, 2), B, c, x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(Eigen::MatrixXd::Ones(2, 5), B, c, x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(A, Eigen::MatrixXd::Ones(5, 1), c, x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(A, B, Eigen::VectorXd::Ones(5), x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(A, B, c, x0, -1), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_WEIGTHS, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);
    auto cost = std::make_shared<copra::TrajectoryCost>(M, xd);

    BOOST_REQUIRE_NO_THROW(cost->weight(2));
    BOOST_REQUIRE_THROW(cost->weights(Eigen::VectorXd::Ones(5)), std::domain_error);
    BOOST_REQUIRE_NO_THROW(cost->weights(wx));
    BOOST_REQUIRE_NO_THROW(controller.addCost(cost));
    BOOST_REQUIRE_NO_THROW(cost->weights(Eigen::VectorXd::Ones(2)));
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_TRAJECTORY_COST, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::TrajectoryCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::TrajectoryCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addCost(badCost2), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_TARGET_COST, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::TargetCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::TargetCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addCost(badCost2), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_CONTROL_COST, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::ControlCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::ControlCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addCost(badCost2), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_MIXED_COST, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badCost1 = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(2, 1), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addCost(badCost1), std::domain_error);
    auto badCost2 = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Identity(2, 1), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addCost(badCost2), std::domain_error);
    auto badCost3 = std::make_shared<copra::MixedCost>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addCost(badCost3), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_TRAJECTORY_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr = std::make_shared<copra::TrajectoryConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr), std::domain_error);
    auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addConstraint(trajConstr), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_CONTROL_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr1 = std::make_shared<copra::ControlConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<copra::ControlConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr2), std::domain_error);

    auto goodConstr = std::make_shared<copra::ControlConstraint>(G, h);
    controller.addConstraint(goodConstr);
    BOOST_REQUIRE_THROW(controller.addConstraint(goodConstr), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_MIXED_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr1 = std::make_shared<copra::MixedConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(2, 1), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<copra::MixedConstraint>(Eigen::MatrixXd::Identity(2, 1), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr2), std::domain_error);
    auto badConstr3 = std::make_shared<copra::MixedConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr3), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_TRAJECTORY_BOUND_CONSTRAINT, BoundedSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr = std::make_shared<copra::TrajectoryBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr), std::domain_error);
    auto tbConstr = std::make_shared<copra::TrajectoryBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(3));
    BOOST_REQUIRE_THROW(controller.addConstraint(tbConstr), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_CONTROL_BOUND_CONSTRAINT, BoundedSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    auto badConstr1 = std::make_shared<copra::ControlBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<copra::ControlBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(3));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr2), std::domain_error);

    auto goodConstr = std::make_shared<copra::ControlBoundConstraint>(uLower, uUpper);
    controller.addConstraint(goodConstr);
    BOOST_REQUIRE_THROW(controller.addConstraint(goodConstr), std::runtime_error);
}

/********************************************************************************************************
 *                               Check remove functions                                                 *
 ********************************************************************************************************/

BOOST_FIXTURE_TEST_CASE(REMOVE_COST_AND_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<copra::PreviewSystem>(A, B, c, x0, nbStep);
    auto controller = copra::LMPC(ps);

    {
        auto xCost = std::make_shared<copra::TargetCost>(M, xd);
        auto uCost = std::make_shared<copra::ControlCost>(N, ud);
        auto trajConstr = std::make_shared<copra::TrajectoryConstraint>(E, f);
        auto contConstr = std::make_shared<copra::ControlConstraint>(G, h);

        controller.addCost(xCost);
        controller.addCost(uCost);
        controller.addConstraint(trajConstr);
        controller.addConstraint(contConstr);

        controller.removeCost(xCost);
        controller.removeCost(uCost);
        controller.removeConstraint(trajConstr);
        controller.removeConstraint(contConstr);
    }

    BOOST_TEST_MESSAGE("\nIn DEBUG mode, if a message appears between this line\n*******");
    controller.solve();
    BOOST_TEST_MESSAGE("*******\nand this line, the remove methods have failed!");
}
