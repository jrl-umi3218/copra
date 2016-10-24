//This file is part of ModelPreviewController.

//ModelPreviewController is free software: you can redistribute it and/or modify
//it under the terms of the GNU Lesser General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//ModelPreviewController is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU Lesser General Public License for more details.
//
//You should have received a copy of the GNU Lesser General Public License
//along with ModelPreviewController.  If not, see <http://www.gnu.org/licenses/>.

#include "SolverInterface.h"

// std
#include <iostream>
#include <utility>

namespace mpc
{

/*
 * SolverInterface
 */
int SolverInterface::SI_fail() const
{
    return 0;
}

int SolverInterface::SI_iter() const
{
    return 0;
}

void SolverInterface::SI_inform() const
{
    std::cout << "No inform() function for this qp" << std::endl;
}

void SolverInterface::SI_printLevel(int /* pl */) const
{
    std::cout << "No printLevel() function for this qp" << std::endl;
}

void SolverInterface::SI_tol(double /* tol */) const
{
    std::cout << "No tol(double) function for this qp" << std::endl;
}

bool SolverInterface::SI_warmStart() const
{
    std::cout << "No warmStart() function for this qp" << std::endl;
    return false;
}

void SolverInterface::SI_warmStart(bool /* w */)
{
    std::cout << "No warmStart(bool) function for this qp" << std::endl;
}

const Eigen::VectorXd &SolverInterface::SI_result() const
{
    return std::move(Eigen::VectorXd());
}

void SolverInterface::SI_problem(int /* nrVar */, int /* nrEq */,
                                 int /* nrInEq */)
{
}

bool SolverInterface::SI_solve(
    const Eigen::MatrixXd & /* Q */, const Eigen::VectorXd & /* C */,
    const Eigen::MatrixXd & /* Aeq */, const Eigen::VectorXd & /* Beq */,
    const Eigen::MatrixXd & /* Aineq */, const Eigen::VectorXd & /* Bineq */,
    const Eigen::VectorXd & /* XL */, const Eigen::VectorXd & /* XU */)
{
    return false;
}

} // namespace pc