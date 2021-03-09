/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "SolverInterface.h"
#include <iostream>
#include <utility>

namespace copra {

/*
 * SolverInterface
 */

int SolverInterface::SI_iter() const
{
    return 0;
}

int SolverInterface::SI_maxIter() const
{
    std::cout << "No maxIter() function for this qp" << std::endl;
    return -1;
}

void SolverInterface::SI_maxIter(int /* maxIter */)
{
    std::cout << "No maxIter(int) function for this qp" << std::endl;
}

void SolverInterface::SI_printLevel(int /* pl */)
{
    std::cout << "No printLevel(int) function for this qp" << std::endl;
}

double SolverInterface::SI_feasibilityTolerance() const
{
    std::cout << "No feasibilityTolerance() function for this qp" << std::endl;
    return -1.;
}

void SolverInterface::SI_feasibilityTolerance(double /* tol */)
{
    std::cout << "No feasibilityTolerance(double) function for this qp" << std::endl;
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

} // namespace copra
