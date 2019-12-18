/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <Eigen/Core>

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

// The final point of the trajectory should be [val, 0] where val can be any value inferior to 0;
// Test bound constraints
struct BoundedSystem {
    BoundedSystem()
        : T(0.005)
        , mass(5)
        , nbStep(300)
        , A(2, 2)
        , B(2, 1)
        , M(2, 2)
        , N(1, 1)
        , c(2)
        , uLower(1)
        , uUpper(1)
        , xLower(2)
        , xUpper(2)
        , x0(2)
        , xd(2)
        , ud(1)
        , wx(2)
        , wu(1)
    {
        // System
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << (-9.81 / 2.) * T * T, -9.81 * T;
        x0 << 0, -5;
        wx << 10, 10000;
        wu << 1e-4;

        // Cost
        M << 1, 0, 0, 1;
        N << 1;
        xd << 0, -1;
        ud << 2;

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
    Eigen::MatrixXd A, B, M, N;
    Eigen::VectorXd c, uLower, uUpper, xLower, xUpper, x0, xd, ud, wx, wu;
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
        , M(2, 2)
        , N(1, 1)
        , c(2)
        , h(1)
        , f(1)
        , x0(2)
        , xd(2)
        , ud(1)
        , wx(2)
        , wu(1)
    {
        // System
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << (-9.81 / 2.) * T * T, -9.81 * T;
        G << 1;
        h << 200; // The force can't be superior to 200
        E << 0, 1;
        f << 0; // The velocity can't be positive
        x0 << 0, -5;
        wx << 10, 10000;
        wu << 1e-4;

        // Cost
        M << 1, 0, 0, 1;
        N << 1;
        xd << 0, -1;
        ud << 2;
    }

    double T, mass;
    int nbStep;
    Eigen::MatrixXd A, B, G, E, M, N;
    Eigen::VectorXd c, h, f, x0, xd, ud, wx, wu;
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
        , M(2, 2)
        , N(1, 1)
        , c(2)
        , f(1)
        , x0(2)
        , xd(2)
        , ud(1)
        , wx(2)
        , wu(1)
    {
        // System
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << (-9.81 / 2.) * T * T, -9.81 * T;
        G << 1;
        E << 0, 1;
        f << 200;
        x0 << 0, -5;
        wx << 10, 10000;
        wu << 1e-4;

        // Cost
        M << 1, 0, 0, 1;
        N << 1;
        xd << 0, -1;
        ud << 2;
    }

    double T, mass;
    int nbStep;
    Eigen::MatrixXd A, B, G, E, M, N;
    Eigen::VectorXd c, f, x0, xd, ud, wx, wu;
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
        , M(2, 2)
        , N(1, 1)
        , c(2)
        , f(2)
        , x0(2)
        , xd(2)
        , ud(1)
        , wx(2)
        , wu(1)
    {
        // System
        A << 1, T, 0, 1;
        B << 0.5 * T * T / mass, T / mass;
        c << (-9.81 / 2.) * T * T, -9.81 * T;
        x0 << 0, 0;
        wx << 10, 10000;
        wu << 1e-4;

        // Cost
        M << 1, 0, 0, 1;
        N << 1;
        xd << 0, 0;
        ud << 2;

        // Trajectory equality constraint
        E.setZero();
        E(0, 0) = 1;
        f = x0;
    }

    double T, mass;
    int nbStep;
    Eigen::MatrixXd A, B, E, M, N;
    Eigen::VectorXd c, f, x0, xd, ud, wx, wu;
};
