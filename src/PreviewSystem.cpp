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

// header
#include "PreviewSystem.h"

//mpc
#include "debugUtils.h"

namespace mpc {

void PreviewSystem::system(const Eigen::MatrixXd& state, const Eigen::MatrixXd& control,
    const Eigen::VectorXd& xInit, const Eigen::VectorXd& xTraj, int numberOfSteps)
{
    system(state, control, Eigen::VectorXd::Zero(state.rows()), xInit, xTraj, numberOfSteps);
}

void PreviewSystem::system(const Eigen::MatrixXd& state, const Eigen::MatrixXd& control,
    const Eigen::VectorXd& bias, const Eigen::VectorXd& xInit,
    const Eigen::VectorXd& xTraj, int numberOfSteps)
{
    if (xInit.rows() != state.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "state", xInit, state));
    if (state.rows() != state.cols())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnSquareMat("state", state));
    if (xInit.rows() != control.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "control", xInit, control));
    if (xInit.rows() != bias.rows())
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRows("xInit", "state", xInit, bias));
    if (numberOfSteps <= 0) // This should desappear later
        DOMAIN_ERROR_EXCEPTION("The number of step sould be a positive number! ");

    isUpdated = false;
    nrUStep = numberOfSteps;
    nrXStep = numberOfSteps + 1;
    xDim = static_cast<int>(state.cols());
    uDim = static_cast<int>(control.cols());
    fullXDim = xDim * nrXStep;
    fullUDim = uDim * nrUStep;

    if (xTraj.rows() != xDim && xTraj.rows() != fullXDim)
        DOMAIN_ERROR_EXCEPTION(throwMsgOnRowsOnPSXDim("xTraj", xTraj, this));

    x0 = xInit;
    xd.resize(fullXDim);
    A = state;
    B = control;
    d = bias;
    Phi.resize(fullXDim, xDim);
    Psi.resize(fullXDim, fullUDim);
    xi.resize(fullXDim);

    auto xTrajLen = static_cast<int>(xTraj.rows());
    for (auto i = 0; i < fullXDim; i += xTrajLen)
        xd.segment(i, xTrajLen) = xTraj;
    Phi.setZero();
    Phi.block(0, 0, xDim, xDim) = Eigen::MatrixXd::Identity(xDim, xDim);
    Psi.setZero();
    xi.setZero();
}

} // namespace mpc