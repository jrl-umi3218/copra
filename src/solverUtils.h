#include "Solver.h"

namespace pc
{

enum class SolverFlag
{
    DEFAULT,
    LSSOL,
    QLD,
    QuadProgDense,
    QuadProgSparse
};

enum class PCFlag
{
    Full,
    Last
};
    
static const std::unique_ptr<SolverInterface> solverFactory(SolverFlag)
{
    switch(flag)
    {
        case SolverFlag::DEFAULT:
            return std::make_unique<QuadProgDenseSolver>();
        case SolverFlag::LSSOL:
            return std::make_unique<LSSOLSolver>();
        case SolverFlag::QLD:
            return std::make_unique<QLDSolver>();
        case SolverFlag::QuadProgDense:
            return std::make_unique<QuadProgDenseSolver>();
        case SolverFlag::QuadProgSparse:
            return std::make_unique<QuadProgSparseSolver>();
    }
}

}