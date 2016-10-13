#include "previewController.h"

namespace pc {

PreviewController::PreviewController(SolverFlag flag = SolverFlag::DEFAULT)
    sol_(pb::solverFactory(flag)),

{
}

} // namespace pc