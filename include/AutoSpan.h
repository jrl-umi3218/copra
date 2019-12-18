/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "typedefs.h"
#include <Eigen/Core>

namespace copra {

/**
 * \brief A class made of static helper function.
 * It can't be instantiated. 
 * This class has helper functions for automatically extending a matrix to a given dimension.
 * The result is a block diagonal matrix.
 */
struct COPRA_DLLAPI AutoSpan {
public:
    // Delete the default constructor. This class should not be instantiated
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
