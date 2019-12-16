// This file is part of copra.

// copra is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// copra is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with copra.  If not, see
// <http://www.gnu.org/licenses/>.

#include "solverUtils.h"

namespace copra {

std::unique_ptr<SolverInterface> solverFactory(SolverFlag flag)
{
    switch (flag) {
#ifdef EIGEN_LSSOL_FOUND
    case SolverFlag::LSSOL:
        return std::unique_ptr<LSSOLSolver>(new LSSOLSolver());
#endif
#ifdef EIGEN_GUROBI_FOUND
    case SolverFlag::GUROBIDense:
        return std::unique_ptr<GUROBISolver>(new GUROBISolver());
#endif
#ifdef EIGEN_QLD_FOUND
    case SolverFlag::QLD:
        return std::unique_ptr<QLDSolver>(new QLDSolver());
#endif
#ifdef EIGEN_OSQP_FOUND
    case SolverFlag::OSQP:
        return std::unique_ptr<OSQPSolver>(new OSQPSolver());
#endif
    // case SolverFlag::QuadProgSparse:
    //    return std::make_unique<QuadProgSparseSolver>();
    case SolverFlag::QuadProgDense:
    default:
        return std::unique_ptr<QuadProgDenseSolver>(new QuadProgDenseSolver());
    }
}

SolverInterface* pythonSolverFactory(SolverFlag flag)
{
    switch (flag) {
#ifdef EIGEN_LSSOL_FOUND
    case SolverFlag::LSSOL:
        return new LSSOLSolver;
#endif
#ifdef EIGEN_GUROBI_FOUND
    case SolverFlag::GUROBIDense:
        return new GUROBISolver;
#endif
#ifdef EIGEN_QLD_FOUND
    case SolverFlag::QLD:
        return new QLDSolver;
#endif
#ifdef EIGEN_OSQP_FOUND
    case SolverFlag::OSQP:
        return new OSQPSolver();
#endif
    case SolverFlag::QuadProgDense:
    default:
        return new QuadProgDenseSolver;
    }
}

} // namespace copra