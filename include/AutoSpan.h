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

#pragma once

// Eigen
#include <Eigen/Core>

// copra
#include "typedefs.h"

namespace copra {

/**
 * \brief A class made of static helper function.
 * It can't be instantiated. 
 * This class has helper functions for automatically extending a matrix to a given dimension.
 * The result is a block diagonal matrix.
 */
class AutoSpan
{
public:
    // Delete the default constrcutor. This class should not be instantiated
    AutoSpan() = delete;

    /**
     * \brief Repeat a matrix until it has new_dim rows.
     * \param mat The matrix to extend. The result is stock in this parameter.
     * \param new_dim The new dimension of the matrix. Must be a multiple of mat.rows().
     * \param addCols Additional zero columns to add at the end of the matrix: the number of columns added is addCols * mat.cols().
     * \throw std::domain_error if new_dim is not a multiple of mat.rows()
     * \note Does nothing if mat.rows() == new_dim
     */ 
    static void spanMatrix(Eigen::MatrixXd& mat, Eigen::Index new_dim, int addCols = 0);

    /**
     * \brief Repeat a vector until its size is new_dim.
     * \param vec The vector to extend. The result is stock in this parameter.
     * \param new_dim The new dimension of the vector. Must be a multiple of vec.size().
     * \throw std::domain_error if new_dim is not a multiple of vec.rows()
     * \note Does nothing if mat.rows() == new_dim
     */ 
    static void spanVector(Eigen::VectorXd& vec, Eigen::Index new_dim);
};

} // namespace copra
