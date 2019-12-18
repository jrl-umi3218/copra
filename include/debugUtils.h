/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include "api.h"

#include "PreviewSystem.h"
#include "typedefs.h"
#include <Eigen/Core>
#include <iostream>

#ifdef NDEBUG
#define CONSTRAINT_DELETION_WARN(warn, format, ...) (void)warn
#else
#define CONSTRAINT_DELETION_WARN(warn, format, ...) \
    if (warn) fprintf(stderr, format, __VA_ARGS__)
#endif

namespace copra {

//see http://stackoverflow.com/questions/37181621/easy-way-of-constructing-information-message-for-throwing-stdexception-using-p
template <class E>
[[noreturn]] inline void COPRA_DLLAPI fancy_throw(std::string msg, char const* file, char const* function, std::size_t line)
{
    throw E(std::string("In file: ") + file + "(line " + std::to_string(line) + "): [In function: " + function + "]\n" + msg);
}

} // namespace copra

#define EXCEPTION(TYPE, MESSAGE) fancy_throw<TYPE>(MESSAGE, __FILE__, __func__, __LINE__)
#define DOMAIN_ERROR_EXCEPTION(MESSAGE) EXCEPTION(std::domain_error, MESSAGE)
#define RUNTIME_ERROR_EXCEPTION(MESSAGE) EXCEPTION(std::runtime_error, MESSAGE)

namespace copra {

/**
 * \brief Message error where two matrices have the same number of rows.
 * \param mat1Name Name of the first checked matrix
 * \param mat2Name Name of the second checked matrix
 * \param mat1 The first checked matrix
 * \param mat2 The second checked matrix
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnRows(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);
/**
 * \brief Message error where two matrices have the same number of rows.
 * This is the same as \see throwMsgOnRows with a longer message.
 * \param mat1Name Name of the first checked matrix
 * \param mat2Name Name of the second checked matrix
 * \param mat1 The first checked matrix
 * \param mat2 The second checked matrix
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnRowsAskAutoSpan(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2);
/**
 * \brief Message error where a matrix is square.
 * \param matName Name of the checked matrix
 * \param mat checked matrix
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnSquareMat(const char* matName, const Eigen::MatrixXd& mat);
/**
 * \brief Message error where a matrix has a number of rows equal to dim
 * \param matName Name of the checked matrix
 * \param mat checked matrix
 * \param dim Value to be equal to
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnRowsOnDim(const char* matName, const Eigen::MatrixXd& mat, Eigen::Index dim);
/**
 * \brief Message error where a matrix has its number of rows equal to preview system xDim
 * \param matName Name of the checked matrix
 * \param mat checked matrix
 * \param ps The preview system
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnRowsOnPSxDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Message error where a matrix has its number of rows equal to preview system xDim or fullXDim
 * \param matName Name of the checked matrix
 * \param mat checked matrix
 * \param ps The preview system
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnRowsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Message error where a matrix has its number of rows equal to preview system uDim or fullUDim
 * \param matName Name of the checked matrix
 * \param mat checked matrix
 * \param ps The preview system
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnRowsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Message error where a matrix has its number of columns equal to preview system xDim or fullXDim
 * \param matName Name of the checked matrix
 * \param mat checked matrix
 * \param ps The preview system
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnColsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Message error where a matrix has its number of columns equal to preview system uDim or fullUDim
 * \param matName Name of the checked matrix
 * \param mat checked matrix
 * \param ps The preview system
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnColsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps);
/**
 * \brief Message error where two matrices number of columns are equal to the preview system datas.
 * Message error where the first matrix and the second matrix are respectively equal to the preview system xDim and uDim.
 * \param mat1Name Name of the first checked matrix
 * \param mat2Name Name of the second checked matrix
 * \param mat1 The first checked matrix
 * \param mat2 The second checked matrix
 * \param ps The preview system
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnColsOnPSxuDim(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2, const PreviewSystem* ps);
/**
 * \brief Message error where two matrices number of columns are equal to the preview system datas.
 * Message error where the first matrix and the second matrix are either respectively equal to the preview system xDim and uDim
 * either respectively equal to the preview system fullXDim and fullUDim.
 * \param mat1Name Name of the first checked matrix
 * \param mat2Name Name of the second checked matrix
 * \param mat1 The first checked matrix
 * \param mat2 The second checked matrix
 * \param ps The preview system
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnColsOnPSXUDim(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2, const PreviewSystem* ps);

/**
 * \brief Message error where a matrix can not be extended to fit the new dimension.
 * \param mat The checked matrix
 * \param new_dim The new dimension asked
 * \return The throwing message
 */
COPRA_DLLAPI std::string throwMsgOnBadNewDim(const Eigen::MatrixXd& mat, Eigen::Index new_dim);
}
