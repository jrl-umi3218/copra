#pragma once

#include "solverUtils.h"
#include <Eigen/Core>
#include <boost/timer/timer.hpp>
#include <memory>
#include <vector>

namespace pc
{

//Forward declaration
class Constrain;

struct PreviewSystemData
{
    PreviewSystemData(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                      const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                      int numberOfSteps);
    int nrStep;
    int xDim, uDim, fullXDim, fullUDim;
    Eigen::VectorXd x0;
    Eigen::VectorXd xd;
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::VectorXd d;
    Eigen::MatrixXd Phi;
    Eigen::MatrixXd Psi;
    Eigen::VectorXd xi;
};

class PreviewController
{
  public:
    PreviewController(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                      const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                      int numberOfSteps, PCFlag pcFlag = PCFlag::Last, SolverFlag sFlag = SolverFlag::DEFAULT);
    PreviewController(const Eigen::MatrixXd &state, const Eigen::MatrixXd &control,
                      const Eigen::VectorXd &bias, const Eigen::VectorXd &xInit, const Eigen::VectorXd &xTraj,
                      int numberOfSteps, PCFlag pcFlag = PCFlag::Last, SolverFlag sFlag = SolverFlag::DEFAULT);

    void selectQPSolveur(SolverFlag flag);
    bool solve();

    const Eigen::VectorXd &control() const noexcept;
    Eigen::VectorXd trajectory() const noexcept;
    boost::timer::cpu_times solveTime() const noexcept;
    boost::timer::cpu_times solveAndBuildTime() const noexcept;

    void weights(double wx, double wu);
    void weights(const Eigen::VectorXd &Wx, const Eigen::VectorXd &Wu);

    void addConstrain(Constrain &constr);
    void resetConstrains() noexcept;

  private:
    void previewSystem();
    void makeQPForm();

  private:
    PCFlag flag_;
    int nrConstr_;

    std::unique_ptr<PreviewSystemData> psd_;
    std::vector<Constrain *> constr_;
    std::unique_ptr<SolverInterface> sol_;

    Eigen::MatrixXd Q_, AInEq_;
    Eigen::VectorXd c_, bInEq_, Wx_, Wu_;

    boost::timer::cpu_timer solveTime_, solveAndBuildTime_;
};

class Constrain
{
  public:
    Constrain(const Eigen::MatrixXd &E, const Eigen::VectorXd &f);

    virtual void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) = 0;
    virtual void update(const std::unique_ptr<PreviewSystemData> &ps) = 0;

    int nrConstr() noexcept
    {
        return nrConstr_;
    }
    const Eigen::MatrixXd &A() noexcept
    {
        return A_;
    }
    const Eigen::VectorXd &b() noexcept
    {
        return b_;
    }

  protected:
    int nrConstr_;
    Eigen::MatrixXd E_, A_;
    Eigen::VectorXd f_, b_;
};

class TrajectoryConstrain final : public Constrain
{
  public:
    void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) override;
    void update(const std::unique_ptr<PreviewSystemData> &psd) override;
};

class ControlConstrain final : public Constrain
{
  public:
    void initializeConstrain(const std::unique_ptr<PreviewSystemData> &psd) override;
    void update(const std::unique_ptr<PreviewSystemData> &ps) override;
};

} // namespace pc