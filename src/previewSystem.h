#pragma once

#include <Eigen/Core>
#include <vector>
#include <memory>
#include "solverUtils.h"

namespace pc
{

struct PreviewSystemData
{
    PreviewSystemData(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                      const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                      std::size_t numberOfSteps);
    std::size_t nrStep;
    std::size_t xDim, uDim, fullXDim, fullUDim;
    Eigen::VectorXd x0;
    Eigen::VectorXd xd;
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::VectorXd d;
    Eigen::MatrixXd Phi;
    Eigen::MatrixXd Psi;
    Eigen::VectorXd xi;
};

//TODO: Add weights
class System
{
  public:
    System(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
           const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
           std::size_t numberOfSteps, PCFlag flag = PCFlag::Last);
    System(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
           const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
           std::size_t numberOfSteps, PCFlag flag = PCFlag::Last);

    void previewSystem();
    void makeQPForm();

    void addConstrain(Constrain &contr);
    void resetConstrains() noexecpt;

    const Eigen::MatrixXd &Q() noexecpt
    {
        return Q_;
    }

    const Eigen::VectorXd &C() noexecpt
    {
        return C_;
    }

    const Eigen::MatrixXd &AInEq() noexecpt
    {
        return AInEq_;
    }

    const Eigen::VectorXd &BInEq() noexecpt
    {
        return BInEq_;
    }

  protected:
    PCFlag flag_;

    Eigen::MatrixXd Q_, AInEq_;
    Eigen::VectorXd c_, bInEq_;

    std::unique_ptr<PreviewSystemData> psd_;
    std::vector<Constrain *> constr_;
};

class Constrain
{
  public:
    Constrain(const Eigen::MatrixXd &E, const Eigen::VectorXd f);

    virtual void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) = 0;
    virtual void update(const PCFlag flag,
                        const std::unique_ptr<PreviewSystemData> &ps) = 0;

    const std::size_t nrConstr() noexecpt
    {
        return nrConstr_;
    }
    const Eigen::MatrixXd &A() noexecpt
    {
        return A_;
    }
    const Eigen::MatrixXd &b() noexecpt
    {
        return b_;
    }

  protected:
    std::size_t nrConstr_;
    Eigen::MatrixXd &E_, A_;
    Eigen::VectorXd &f_, b_;
};

class TrajectoryConstrain final : public Constrain
{
  public:
    void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) override;
    void update(const PCFlag flag, const std::unique_ptr<PreviewSystemData> &psd) override;
};

class ControlConstrain final : public Constrain
{
  public:
    void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) override;
    void update(const PCFlag flag, const std::unique_ptr<PreviewSystemData> &ps) override;
};
}