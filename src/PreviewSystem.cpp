/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "PreviewSystem.h"
#include "debugUtils.h"

namespace copra {

PreviewSystem::PreviewSystem(const Eigen::MatrixXd& state, const Eigen::MatrixXd& control,
    const Eigen::VectorXd& bias, const Eigen::VectorXd& xInit, int numberOfSteps)
{
    system(state, control, bias, xInit, numberOfSteps);
}

void PreviewSystem::system(const Eigen::MatrixXd& state, const Eigen::MatrixXd& control,
    const Eigen::VectorXd& bias, const Eigen::VectorXd& xInit, int numberOfSteps)
{
    if (xInit.rows() != state.rows()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "state", xInit, state));
    }
    if (state.rows() != state.cols()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnSquareMat("state", state));
    }
    if (xInit.rows() != control.rows()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "control", xInit, control));
    }
    if (xInit.rows() != bias.rows()) {
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "state", xInit, bias));
    }
    if (numberOfSteps <= 0) {
        DOMAIN_ERROR_EXCEPTION("The number of step sould be a positive number! ");
    }

    isUpdated = false;
    nrUStep = numberOfSteps;
    nrXStep = numberOfSteps + 1;
    xDim = static_cast<int>(state.cols());
    uDim = static_cast<int>(control.cols());
    fullXDim = xDim * nrXStep;
    fullUDim = uDim * nrUStep;

    x0 = xInit;
    A = state;
    B = control;
    d = bias;
    Phi.resize(fullXDim, xDim);
    Psi.resize(fullXDim, fullUDim);
    xi.resize(fullXDim);

    Phi.setZero();
    Phi.block(0, 0, xDim, xDim) = Eigen::MatrixXd::Identity(xDim, xDim);
    Psi.setZero();
    xi.setZero();
}

void PreviewSystem::updateSystem() noexcept
{
    Phi.block(xDim, 0, xDim, xDim) = A;
    Psi.block(xDim, 0, xDim, uDim) = B;
    xi.segment(xDim, xDim) = d;

    for (auto i = 2; i < nrXStep; ++i) {
        Phi.block(i * xDim, 0, xDim, xDim).noalias() = A * Phi.block((i - 1) * xDim, 0, xDim, xDim);
        Psi.block(i * xDim, 0, xDim, uDim).noalias() = A * Psi.block((i - 1) * xDim, 0, xDim, uDim);
        for (auto j = 1; j < i; ++j) {
            Psi.block(i * xDim, j * uDim, xDim, uDim).noalias() = Psi.block((i - 1) * xDim, (j - 1) * uDim, xDim, uDim);
        }

        xi.segment(i * xDim, xDim).noalias() = A * xi.segment((i - 1) * xDim, xDim) + d;
    }

    isUpdated = true;
}

void PreviewSystem::computeTrajectory(const Eigen::VectorXd& control, Eigen::VectorXd& trajectory) noexcept
{
    trajectory.resize(fullXDim);
    trajectory = this->Phi * this->x0 + this->Psi * control + this->xi;
}

} // namespace copra
