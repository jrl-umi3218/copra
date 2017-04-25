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

#define BOOST_TEST_MODULE TestPreviewControl

// stl
#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

// boost
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// eigen
#include <Eigen/Core>

// mpc
#include "MPC.h"
#include "PreviewSystem.h"
#include "constraints.h"
#include "costFunctions.h"
#include "solverConfig.h"
#include "solverUtils.h"

// The final point of the trajectory should be [val, 0] where val can be any value inferior to 0;
// Test bound constraints
struct BoundedSystem {
    BoundedSystem()
        : T(0.005)
        , mass(5)
        , nbStep(300)
        , A(2, 2)
        , B(2, 1)
        , c(2)
        , uLower(1)
        , uUpper(1)
        , xLower(2)
        , xUpper(2)
        , x0(2)
        , xd(2)
        , wx(2)
        , wu(1)
    {
        // System
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << (-9.81 / 2.) * T * T, -9.81 * T;
        x0 << 0, -5;
        xd << 0, 0;
        wx << 10, 10000;
        wu << 1e-4;

        // Control bound
        uLower.setConstant(-std::numeric_limits<double>::infinity());
        uUpper.setConstant(200); // The force can't be superior to 200

        // Trajectory bound
        xLower.setConstant(-std::numeric_limits<double>::infinity());
        xUpper(0) = std::numeric_limits<double>::infinity();
        xUpper(1) = 0; // The velocity can't be positive
    }

    double T, mass;
    int nbStep;
    Eigen::MatrixXd A, B;
    Eigen::VectorXd c, uLower, uUpper, xLower, xUpper, x0, xd, wx, wu;
};

// The final point of the trajectory should be [val, 0] where val can be any value inferior to 0 (same as previous one)
// Test inequality constraints
struct IneqSystem {
    IneqSystem()
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
        c << (-9.81 / 2.) * T * T, -9.81 * T;
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

// The final point of the trajectory should be [val, 0] where val can be any value inferior to 0 (same as previous one)
// Test mixed constraints
struct MixedSystem {
    MixedSystem()
        : T(0.005)
        , mass(5)
        , nbStep(300)
        , A(2, 2)
        , B(2, 1)
        , G(1, 1)
        , E(1, 2)
        , c(2)
        , f(1)
        , x0(2)
        , xd(2)
        , wx(2)
        , wu(1)
    {
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << (-9.81 / 2.) * T * T, -9.81 * T;
        G << 1;
        E << 0, 1;
        f << 200;
        x0 << 0, -5;
        xd << 0, 0;
        wx << 10, 10000;
        wu << 1e-4;
    }

    double T, mass;
    int nbStep;
    Eigen::MatrixXd A, B, G, E;
    Eigen::VectorXd c, f, x0, xd, wx, wu;
};

// Search forces that let the system immobile (should be equal to gravity * tiemstep)
// Test Equality constraints
// xd becomes useless here
struct EqSystem {
    EqSystem()
        : T(0.005)
        , mass(5)
        , nbStep(300)
        , A(2, 2)
        , B(2, 1)
        , E(2, 2)
        , c(2)
        , f(2)
        , x0(2)
        , xd(2)
        , wx(2)
        , wu(1)
    {
        // System
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << (-9.81 / 2.) * T * T, -9.81 * T;
        x0 << 0, 0;
        xd << 0, 0;
        wx << 10, 10000;
        wu << 1e-4;

        // Trajectory equality constraint
        E.setZero();
        E(0, 0) = 1;
        f = x0;
    }

    bool printErrorMessage(const std::domain_error& e)
    {
        std::cout << e.what() << std::endl
                  << std::endl;
        return true;
    }

    double T, mass;
    int nbStep;
    Eigen::MatrixXd A, B, E;
    Eigen::VectorXd c, f, x0, xd, wx, wu, uLower, uUpper;
};

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_BOUND_CONSTRAINTS, BoundedSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TargetCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto trajConstr = std::make_shared<mpc::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<mpc::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";

    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_BOUND_CONSTRAINTS, BoundedSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TrajectoryCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto trajConstr = std::make_shared<mpc::TrajectoryBoundConstraint>(xLower, xUpper);
    auto contConstr = std::make_shared<mpc::ControlBoundConstraint>(uLower, uUpper);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";
    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_BOUND_CONSTRAINT, BoundedSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::VectorXd& xLower, const Eigen::VectorXd& xUpper, const Eigen::VectorXd& uLower, const Eigen::VectorXd& uUpper) {
        auto trajConstr = std::make_shared<mpc::TrajectoryBoundConstraint>(xLower, xUpper);
        trajConstr->autoSpan();

        auto contConstr = std::make_shared<mpc::ControlBoundConstraint>(uLower, uUpper);
        contConstr->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addConstraint(trajConstr));
        BOOST_REQUIRE_NO_THROW(controller.addConstraint(contConstr));
    };

    auto spanVector = [&](const Eigen::VectorXd& v, int size) {
        Eigen::VectorXd vout = Eigen::VectorXd::Zero(v.rows() * size);
        for (int i = 0; i < size; ++i)
            vout.segment(i * v.rows(), v.rows()) = v;
        return vout;
    };

    auto fullxLower = spanVector(xLower, nbXStep);
    auto fullxUpper = spanVector(xUpper, nbXStep);
    auto fulluLower = spanVector(uLower, nbStep);
    auto fulluUpper = spanVector(uUpper, nbStep);

    checkSpan(fullxLower, xUpper, fulluLower, uUpper);
    checkSpan(xLower, fullxUpper, uLower, fulluUpper);
    checkSpan(fullxLower, fullxUpper, fulluLower, fulluUpper);
}

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_INEQUALITY_CONSTRAINTS, IneqSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TargetCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(E, f);
    auto contConstr = std::make_shared<mpc::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";

    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_INEQUALITY_CONSTRAINTS, IneqSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TrajectoryCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(E, f);
    auto contConstr = std::make_shared<mpc::ControlConstraint>(G, h);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);
    controller.addConstraint(contConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";
    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_INEQUALITY_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    int nbXStep = nbStep + 1;

    auto checkSpan = [&](const Eigen::MatrixXd& E, const Eigen::VectorXd& f, const Eigen::MatrixXd& G, const Eigen::VectorXd& h) {
        auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(E, f);
        trajConstr->autoSpan();

        auto contConstr = std::make_shared<mpc::ControlConstraint>(G, h);
        contConstr->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addConstraint(trajConstr));
        BOOST_REQUIRE_NO_THROW(controller.addConstraint(contConstr));
    };

    auto spanMatrix = [&](const Eigen::MatrixXd& m, int size) {
        Eigen::MatrixXd mout = Eigen::MatrixXd::Zero(m.rows() * size, m.cols() * size);
        for (int i = 0; i < size; ++i)
            mout.block(i * m.rows(), i * m.cols(), m.rows(), m.cols()) = m;
        return mout;
    };
    auto spanVector = [&](const Eigen::VectorXd& v, int size) {
        Eigen::VectorXd vout = Eigen::VectorXd::Zero(v.rows() * size);
        for (int i = 0; i < size; ++i)
            vout.segment(i * v.rows(), v.rows()) = v;
        return vout;
    };

    auto fullE = spanMatrix(E, nbXStep);
    auto fullf = spanVector(f, nbXStep);
    auto fullG = spanMatrix(G, nbStep);
    auto fullh = spanVector(h, nbStep);

    checkSpan(fullE, f, fullG, h);
    checkSpan(E, fullf, G, fullh);
    checkSpan(fullE, fullf, fullG, fullh);
}

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_MIXED_CONSTRAINTS, MixedSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TargetCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto mixedConstr = std::make_shared<mpc::MixedConstraint>(E, G, f);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";

    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_MIXED_CONSTRAINTS, MixedSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TrajectoryCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto mixedConstr = std::make_shared<mpc::MixedConstraint>(E, G, f);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(mixedConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";
    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(CHECK_AUTOSPAN_AND_WHOLE_MATRIX_ON_MIXED_CONSTRAINT, MixedSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);

    auto checkSpan = [&](const Eigen::MatrixXd& E, const Eigen::MatrixXd& G, const Eigen::VectorXd& f) {
        auto mixedConstr = std::make_shared<mpc::MixedConstraint>(E, G, f);
        mixedConstr->autoSpan();

        BOOST_REQUIRE_NO_THROW(controller.addConstraint(mixedConstr));
    };

    auto spanMatrix = [&](const Eigen::MatrixXd& m, int size, int addCols = 0) {
        Eigen::MatrixXd mout = Eigen::MatrixXd::Zero(m.rows() * size, m.cols() * (size + addCols));
        for (int i = 0; i < size; ++i)
            mout.block(i * m.rows(), i * m.cols(), m.rows(), m.cols()) = m;
        return mout;
    };
    auto spanVector = [&](const Eigen::VectorXd& v, int size) {
        Eigen::VectorXd vout = Eigen::VectorXd::Zero(v.rows() * size);
        for (int i = 0; i < size; ++i)
            vout.segment(i * v.rows(), v.rows()) = v;
        return vout;
    };

    auto fullE = spanMatrix(E, nbStep, 1);
    auto fullG = spanMatrix(G, nbStep);
    auto fullf = spanVector(f, nbStep);

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

BOOST_FIXTURE_TEST_CASE(MPC_TARGET_COST_WITH_EQUALITY_CONSTRAINTS, EqSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TargetCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(E, f, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";

    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(MPC_TRAJECTORY_COST_WITH_EQUALITY_CONSTRAINTS, EqSystem)
{
    std::vector<std::pair<std::string, double> > solveTime;

    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);
    auto xCost = std::make_shared<mpc::TrajectoryCost>(Eigen::MatrixXd::Identity(2, 2), xd);
    auto uCost = std::make_shared<mpc::ControlCost>(Eigen::MatrixXd::Identity(1, 1), Eigen::VectorXd::Zero(1));
    auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(E, f, false);
    xCost->weights(wx);
    uCost->weights(wu);

    controller.addCost(xCost);
    controller.addCost(uCost);
    controller.addConstraint(trajConstr);

    auto pcCheck = [&](const std::string& solverName, mpc::SolverFlag sFlag) {
        controller.selectQPSolver(sFlag);

        BOOST_REQUIRE(controller.solve());
        solveTime.emplace_back(solverName, controller.solveTime() * 1e3);

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

    pcCheck("Default (QuadProgDense)", mpc::SolverFlag::DEFAULT);
#ifdef EIGEN_LSSOL_FOUND
    pcCheck("LSSOL", mpc::SolverFlag::LSSOL);
#endif
#ifdef EIGEN_QLD_FOUND
    pcCheck("QLD", mpc::SolverFlag::QLD);
#endif
#ifdef EIGEN_GUROBI_FOUND
    pcCheck("GUROBIDense", mpc::SolverFlag::GUROBIDense);
#endif

    std::sort(solveTime.begin(), solveTime.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    std::stringstream ss;
    ss << "Solving fasteness: ";
    for (auto sol : solveTime)
        ss << sol.first << " (" + std::to_string(sol.second) << "ms) > ";
    BOOST_MESSAGE(ss.str());
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_PREVIEW_SYSTEM, IneqSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    BOOST_REQUIRE_THROW(ps->system(Eigen::MatrixXd::Ones(5, 2), B, c, x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(Eigen::MatrixXd::Ones(2, 5), B, c, x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(A, Eigen::MatrixXd::Ones(5, 1), c, x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(A, B, Eigen::VectorXd::Ones(5), x0, nbStep), std::domain_error);
    BOOST_REQUIRE_THROW(ps->system(A, B, c, x0, -1), std::domain_error);
}

/*
BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_WEIGTHS, IneqSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto fullController = mpc::MPCTypeFull(ps);
    auto lastController = mpc::MPCTypeLast(ps);

    fullController.weights(1.1, 2.2);
    lastController.weights(1.1, 2.2);
    BOOST_REQUIRE_THROW(fullController.weights(Eigen::VectorXd::Ones(5), wu), std::domain_error);
    BOOST_REQUIRE_THROW(fullController.weights(wx, Eigen::VectorXd::Ones(5)), std::domain_error);
    BOOST_REQUIRE_THROW(lastController.weights(Eigen::VectorXd::Ones(5), wu), std::domain_error);
    BOOST_REQUIRE_THROW(lastController.weights(wx, Eigen::VectorXd::Ones(5)), std::domain_error);
    try {
        lastController.weights(wx, Eigen::VectorXd::Ones(5));
    } catch (const std::domain_error& e) {
        std::cerr << "Test error message output" << std::endl;
        std::cerr << e.what() << std::endl
                  << std::endl;
    }
}
*/

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_TRAJECTORY_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);

    auto badConstr = std::make_shared<mpc::TrajectoryConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr), std::domain_error);
    auto trajConstr = std::make_shared<mpc::TrajectoryConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addConstraint(trajConstr), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_CONTROL_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);

    auto badConstr1 = std::make_shared<mpc::ControlConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<mpc::ControlConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr2), std::domain_error);

    auto goodConstr = std::make_shared<mpc::ControlConstraint>(G, h);
    controller.addConstraint(goodConstr);
    BOOST_REQUIRE_THROW(controller.addConstraint(goodConstr), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_MIXED_CONSTRAINT, IneqSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);

    auto badConstr1 = std::make_shared<mpc::MixedConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(2, 1), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<mpc::MixedConstraint>(Eigen::MatrixXd::Identity(2, 1), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr2), std::domain_error);
    auto badConstr3 = std::make_shared<mpc::MixedConstraint>(Eigen::MatrixXd::Identity(5, 5), Eigen::MatrixXd::Identity(5, 5), Eigen::VectorXd::Ones(5));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr3), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_TRAJECTORY_BOUND_CONSTRAINT, BoundedSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);

    auto badConstr = std::make_shared<mpc::TrajectoryBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr), std::domain_error);
    auto tbConstr = std::make_shared<mpc::TrajectoryBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(3));
    BOOST_REQUIRE_THROW(controller.addConstraint(tbConstr), std::domain_error);
}

BOOST_FIXTURE_TEST_CASE(ERROR_HANDLER_FOR_CONTROL_BOUND_CONSTRAINT, BoundedSystem)
{
    auto ps = std::make_shared<mpc::PreviewSystem>();
    ps->system(A, B, c, x0, nbStep);
    auto controller = mpc::MPC(ps);

    auto badConstr1 = std::make_shared<mpc::ControlBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(2));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr1), std::domain_error);
    auto badConstr2 = std::make_shared<mpc::ControlBoundConstraint>(Eigen::VectorXd::Ones(3), Eigen::VectorXd::Ones(3));
    BOOST_REQUIRE_THROW(controller.addConstraint(badConstr2), std::domain_error);

    auto goodConstr = std::make_shared<mpc::ControlBoundConstraint>(uLower, uUpper);
    controller.addConstraint(goodConstr);
    BOOST_REQUIRE_THROW(controller.addConstraint(goodConstr), std::runtime_error);
}