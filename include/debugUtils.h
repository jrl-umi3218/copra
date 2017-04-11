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

template <class E> //see http://stackoverflow.com/questions/37181621/easy-way-of-constructing-information-message-for-throwing-stdexception-using-p
[[noreturn]] void fancy_throw(std::string msg, char const* file, char const* function, std::size_t line)
{
    throw E(std::string("In file: ") + file + "(line " + std::to_string(line) + "): [In function: " + function + "]\n" + msg);
}
} // namespace mpc

#define EXCEPTION(TYPE, MESSAGE) \
    fancy_throw<TYPE>(MESSAGE, __FILE__, __func__, __LINE__)
#define DOMAIN_ERROR_EXCEPTION(MESSAGE) EXCEPTION(std::domain_error, MESSAGE)

namespace mpc {

/**
 * \brief Check that a matrix has the same dimension of another one.
 * \param isMatName Name of the matrix to check
 * \param isMat Matrix to check
 * \param shouldBeMat Matrix to be compared with
 * \return The throwing message
 */
std::string throwMsgOnMat(const char* isMatName, const Eigen::MatrixXd& isMat, const Eigen::MatrixXd& shouldBeMat);
/**
 * \brief Check that two matrices have the same number of rows.
 * \param mat1Name Name of the first matrix to check
 * \param mat2Name Name of the second matrix to check
 * \param mat1 The first matrix to check
 * \param mat2 The second matrix to check
 * \return The throwing message
 */
std::string throwMsgOnRows(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);
/**
 * \brief Check that a matrix is square.
 * \param matName Name of the matrix to check
 * \param mat Matrix to check
 * \return The throwing message
 */
std::string throwMsgOnSquareMat(const char* matName, const Eigen::MatrixXd& mat);
/**
 * \brief Check that a matrix has a number of rows equal to dim
 * \param matName Name of the matrix to check
 * \param mat Matrix to check
 * \param dim Value to be equal to
 * \return The throwing message
 */
std::string throwMsgOnRowsOnDim(const char* matName, const Eigen::MatrixXd& mat, Eigen::Index dim);
/**
 * \brief Check that a matrix has its number of rows equal to preview system xDim or fullXDim
 * \param matName Name of the matrix to check
 * \param mat Matrix to check
 * \param ps The preview system
 * \return The throwing message
 */
std::string throwMsgOnRowsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Check that a matrix has its number of rows equal to preview system uDim or fullUDim
 * \param matName Name of the matrix to check
 * \param mat Matrix to check
 * \param ps The preview system
 * \return The throwing message
 */
std::string  throwMsgOnRowsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Check that a matrix has its number of columns equal to preview system xDim or fullXDim
 * \param matName Name of the matrix to check
 * \param mat Matrix to check
 * \param ps The preview system
 * \return The throwing message
 */
std::string throwMsgOnColsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Check that a matrix has its number of columns equal to preview system uDim or fullUDim
 * \param matName Name of the matrix to check
 * \param mat Matrix to check
 * \param ps The preview system
 * \return The throwing message
 */
std::string throwMsgOnColsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Check that two matrices number of columns are equal to the preview system datas.
 * Check that the first matrix and the second matrix are either respectively equal to the preview system xDim and uDim
 * either respectively equal to the preview system fullXDim and fullUDim.
 * \param mat1Name Name of the first matrix to check
 * \param mat2Name Name of the second matrix to check
 * \param mat1 The first matrix to check
 * \param mat2 The second matrix to check
 * \param ps The preview system
 * \return The throwing message
 */
std::string throwMsgOnColsOnPSXUDim(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2, const PreviewSystem* ps);
}