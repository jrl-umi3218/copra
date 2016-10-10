#pragma once

#include <Eigen/Core>

namespace pc
{

class System
{
public:
    System(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
        double loopTime, double window);

    virtual void trajectoryConstrain(const Eigen::MatrixXd& E, const Eigen::VectorXd f);
    void controlConstrain(const Eigen::MatrixXd& E, const Eigen::VectorXd f);
    virtual void previewSystem(PreviewControllerFlag flag) noexecpt;
    virtual void makeQPForm() noexecpt;

    const Eigen::MatrixXd& Q() noexecpt
    {
        return Q_;
    }

    const Eigen::VectorXd& C() noexecpt
    {
        return C_;
    }

    const Eigen::MatrixXd& AInEq() noexecpt
    {
        return AInEq_;
    }

    const Eigen::VectorXd& BInEq() noexecpt
    {
        return BInEq_;
    }

protected:
    double loopTime_, window_;
    std::size_t nbStep_;
    int xDim_, uDim_, fullXDim_, fullUDim_;

    Eigen::MatrixXd A_, B_, phi_, psi_, Q_, AInEq_;
    Eigen::VectorXd C_, BInEq_;
};

class SystemWithConstantBias final : public System
{
public:
    SystemWithConstantBias(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B,
        const Eigen::VectorXd& C);

    void trajectoryConstrain(const Eigen::MatrixXd& E, const Eigen::VectorXd f) override;
    const Eigen::MatrixXd& Xi() noexecpt
    {
        return xi_;
    } 

    void previewSystem(PreviewControllerFlag flag) override noexecpt;
    void makeQPForm() override noexecpt;

protected:
    Eigen::VectorXd xi_, d_;    
};

};