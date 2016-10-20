#include "PreviewSystem.h"

namespace mpc
{

void PreviewSystem::system(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                           const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                           int numberOfSteps)
{

    system(state, control, Eigen::VectorXd::Zero(state.rows()), xInit, xTraj, numberOfSteps);
}

void PreviewSystem::system(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                           const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                           int numberOfSteps)
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