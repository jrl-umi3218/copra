#include "previewController.h"


namespace pc
{

PreviewController::PreviewController(SolverFlag flag=SolverFlag::DEFAULT)
    sol_(pb::solverFactory(flag)),

{
}

PreviewController::computePCMatrices()
{
    int i = 0;
    int j = 0;
    auto lineX = [this, i](int var=0)
    {
        return (i + var)*xDim_;
    }
    auto lineU = [this, j](int var=0)
    {
        return (j + var)*uDim_;
    }

    for(i = 1; i < nbStep_; ++i)
    {
        phi_.block(lineX(), 0, xDim_, xDim_) = A_*phi_.block(lineX(-1), 0, xDim_, xDim_);
        for(j = 0; j < nbStep_; ++j)
            psi_.block(lineX(), lineU(), xDim_, uDim_) = A_*psi_.block(lineX(-1), lineU(), xDim_, uDim_)

        psi_.block(lineX(), lineX(), xDim_, uDim_) = B_;
        xi_.segment(lineX(), xDim_) = A_*xi_.segment(lineX(-1), xDim_) + c_; 
    }

    
}

} // namespace pc