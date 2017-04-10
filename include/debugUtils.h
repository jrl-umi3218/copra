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

#pragma once

// stl
#include <iostream>

// Eigen
#include <Eigen/Core>

// mpc
#include "PreviewSystem.h"

#ifndef _DEBUG
#define CONSTRAINT_DELETION_WARN(warn, format, ...)
#else
#define CONSTRAINT_DELETION_WARN(warn, format, ...) \
    if ((warn)) {                                   \
        fprintf(stderr, format, __VA_ARGS__);       \
    }
#endif

namespace mpc {

void checkMat(const char* isMatName, const Eigen::MatrixXd& isMat, const Eigen::MatrixXd& shouldBeMat);
void checkRows(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);
void checkSquaredMat(const char* matName, const Eigen::MatrixXd& mat);
void checkRowsOnDim(const char* mat1Name, const Eigen::MatrixXd& mat1, Eigen::Index dim);
void checkRowsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
void checkRowsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
void checkColsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
void checkColsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
void checkColsOnPSXUDim(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2, const PreviewSystem* ps);
}