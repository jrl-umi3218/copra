#include "SolverInterface.h"
#include "QuadProgSolver.h"
#include "QLDSolver.h"
#ifdef LSSOL_SOLVER_FOUND
#include "LSSOLSolver.h"
#endif

#include <memory>

namespace pc
{

enum class SolverFlag
{
    DEFAULT,
#ifdef LSSOL_SOLVER_FOUND
    LSSOL,
#endif
    QLD,
    QuadProgDense,
    // QuadProgSparse
};

enum class PCFlag
{
    Full,
    Last
};

SolverInterface solverFactory(SolverFlag flag)
{
    switch (flag)
    {
    case SolverFlag::DEFAULT:
        return QuadProgDenseSolver();
#ifdef LSSOL_SOLVER_FOUND
    case SolverFlag::LSSOL:
        return LSSOLSolver();
#endif
    case SolverFlag::QLD:
        return QLDSolver();
    case SolverFlag::QuadProgDense:
        return QuadProgDenseSolver();
        // case SolverFlag::QuadProgSparse:
        //    return QuadProgSparseSolver();
    }
}

} // namespace pc