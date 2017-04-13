// This file is part of mpc.

// mpc is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// mpc is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with mpc.  If not, see
// <http://www.gnu.org/licenses/>.

#include "solverUtils.h"

namespace mpc {

std::unique_ptr<SolverInterface> MPC_DLLAPI solverFactory(SolverFlag flag)
{
    switch (flag) {
#ifdef EIGEN_LSSOL_FOUND
    case SolverFlag::LSSOL:
        return std::unique_ptr<LSSOLSolver>();
#endif
#ifdef EIGEN_GUROBI_FOUND
    case SolverFlag::GUROBIDense:
        return std::unique_ptr<GUROBISolver>();
#endif
#ifdef EIGEN_QLD_FOUND
    case SolverFlag::QLD:
        return std::unique_ptr<QLDSolver>();
#endif
    // case SolverFlag::QuadProgSparse:
    //    return std::make_unique<QuadProgSparseSolver>();
    case SolverFlag::QuadProgDense:
    default:
        return std::unique_ptr<QuadProgDenseSolver>();
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
    case SolverFlag::QuadProgDense:
    default:
        return new QuadProgDenseSolver;
    }
}

} // namespace mpc