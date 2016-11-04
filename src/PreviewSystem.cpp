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

#include "PreviewSystem.h"

namespace mpc {

void PreviewSystem::system(const Eigen::MatrixXd& state,
    const Eigen::MatrixXd& control,
    const Eigen::VectorXd& xInit,
    const Eigen::VectorXd& xTraj, int numberOfSteps)
{
    system(state, control, Eigen::VectorXd::Zero(state.rows()), xInit, xTraj,
        numberOfSteps);
}

void PreviewSystem::system(const Eigen::MatrixXd& state,
    const Eigen::MatrixXd& control,
    const Eigen::VectorXd& bias, const Eigen::VectorXd& xInit,
    const Eigen::VectorXd& xTraj, int numberOfSteps)
{
    assert(state.rows() == control.rows());
    assert(state.rows() == bias.rows());
    assert(state.rows() == xInit.rows());
    assert(state.rows() == xTraj.rows());
    assert(state.cols() == xInit.rows());

    nrStep = numberOfSteps;
    xDim = static_cast<int>(state.cols());
    uDim = static_cast<int>(control.cols());
    fullXDim = xDim * nrStep;
    fullUDim = uDim * nrStep;

    assert(xTraj.rows() == xDim || xTraj.rows() == fullXDim);

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
    Psi.setZero();
    xi.setZero();
}

} // namespace mpc