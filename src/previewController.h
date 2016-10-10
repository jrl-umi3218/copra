#pragma once

// std
#include <memory>

// Eigen
#include <Eigen/Core>

// Solveurs
#include "solverUtils.h"

namespace pc
{

class PreviewController
{
public:


public:
    PreviewController(PreviewControllerFlag flag, SolverFlag flag=SolverFlag::DEFAULT);

    void selectQPSolveur(SolverFlag flag);
    bool solve();
    const Eigen::VectorXd& trajectory(); 

private:
    void computePCMatrices();

private:
    std::unique_ptr<SolverInterface> sol_;
    TypePCFlag flag_;

    std::shared_ptr<System> system_;
};

} // namespace pc