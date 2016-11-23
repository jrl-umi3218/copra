// This file is part of ModelPreviewController.

// ModelPreviewController is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ModelPreviewController is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ModelPreviewController.  If not, see
// <http://www.gnu.org/licenses/>.

#define BOOST_TEST_MODULE TestPreviewControl
#include <Eigen/Core>
#include <algorithm>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <memory>
#include <numeric>
#include <vector>

#include "Constraints.h"
#include "PreviewController.h"
#include "PreviewSystem.h"
#include "solverUtils.h"

// The final point of the trajectory should be [val, 0] where val can be any
// value inferior to 0;
struct System {
    System()
        : T(0.005)
        , mass(5)
        , nbStep(300)
        , A(2, 2)
        , B(2, 1)
        , G(1, 1)
        , E(1, 2)
        , c(2)
        , h(1)
        , f(1)
        , x0(2)
        , xd(2)
        , wx(2)
        , wu(1)
    {
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << 0, -9.81 * T;
        G << 1;
        h << 200; // The force can't be superior to 200
        E << 0, 1;
        f << 0; // The velocity can't be positive
        x0 << 0, -5;
        xd << 0, 0;
        wx << 10, 10000;
        wu << 1e-4;
    }

    double T, mass;
    int nbStep;
    Eigen::MatrixXd A, B, G, E;
    Eigen::VectorXd c, h, f, x0, xd, wx, wu;
};

BOOST_FIXTURE_TEST_CASE(OneDofSystemTypeLast, System)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, xd, nbStep);
    auto controller = mpc::MPCTypeLast(ps);
    auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(E, f);
    auto contConstr = std::make_shared<mpc::ControlConstraint>(G, h);

    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    controller.weights(wx, wu);

    auto pcCheck = [&](const std::string &solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(
            solverName, static_cast<double>(controller.solveTime().wall) * 1e-6);

        Eigen::VectorXd fullTraj = controller.trajectory();
        auto trajLen = fullTraj.rows() / 2;
        Eigen::VectorXd posTraj(trajLen);
        Eigen::VectorXd velTraj(trajLen);
        for (auto i = 0; i < trajLen; ++i)
        {
            posTraj(i) = fullTraj(2 * i);
            velTraj(i) = fullTraj(2 * i + 1);
        }
        Eigen::VectorXd control = controller.control();

        // Check result
        BOOST_CHECK_SMALL(xd(1) - velTraj.tail(1)(0), 0.001);

        // Check constrains
        BOOST_REQUIRE_LE(posTraj.maxCoeff(), x0(0));
        BOOST_REQUIRE_LE(control.maxCoeff(), h(0));
    };

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
    pcCheck("QLD", mpc::SolverFlag::QLD);
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);

    std::sort(
        solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime) {
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";
    }
    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(OneDofSystemTypeFull, System)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, xd, nbStep);
    auto controller = mpc::MPCTypeFull(ps);
    auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(E, f);
    auto contConstr = std::make_shared<mpc::ControlConstraint>(G, h);

    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    controller.weights(wx, wu);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(
            solverName, static_cast<double>(controller.solveTime().wall) * 1e-6);

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
        BOOST_REQUIRE_LE(
            control.maxCoeff(),
            h(0) + 1e-3); // QuadProg allows to exceeds the constrain of a small mount.
    };

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
    pcCheck("QLD", mpc::SolverFlag::QLD);
    pcCheck("GUROBI", mpc::SolverFlag::GUROBIDense);

    std::sort(
        solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime) {
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";
    }
    BOOST_MESSAGE(ss.str());
}