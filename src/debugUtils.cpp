/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "debugUtils.h"
#include <sstream>

namespace copra {

std::string throwMsgOnRows(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2)
{
    std::ostringstream os;
    os << mat1Name << " and " << mat2Name << " should have same number of rows. "
       << mat1Name << " has " << std::to_string(mat1.rows()) << " rows and "
       << mat2Name << " has " << std::to_string(mat2.rows()) << " rows.";
    return os.str();
}

std::string throwMsgOnRowsAskAutoSpan(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2)
{
    std::ostringstream os(throwMsgOnRows(mat1Name, mat2Name, mat1, mat2));
    os << " Maybe you have forgottent to call autoSpanMatrix method?\n";
    return os.str();
}

std::string throwMsgOnSquareMat(const char* matName, const Eigen::MatrixXd& mat)
{
    std::ostringstream os;
    os << matName << " should be a squared matrix but it is a ("
       << mat.rows() << "-by-"
       << mat.cols() << ") matrix.";
    return os.str();
}

std::string throwMsgOnRowsOnDim(const char* matName, const Eigen::MatrixXd& mat, Eigen::Index dim)
{
    std::ostringstream os;
    os << matName << " should have exactly " << std::to_string(dim) << " rows. But "
       << matName << " has " << std::to_string(mat.rows()) << " rows.";
    return os.str();
}

std::string throwMsgOnRowsOnPSxDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    std::ostringstream os;
    os << matName << " should have " << std::to_string(ps->xDim) << " rows. But "
       << matName << " has " << std::to_string(mat.rows()) << " rows.";
    return os.str();
}

std::string throwMsgOnRowsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    std::ostringstream os;
    os << matName << " should have " << std::to_string(ps->xDim) << " or " << std::to_string(ps->fullXDim)
       << " rows. But " << matName << " has " << std::to_string(mat.rows()) << " rows.";
    return os.str();
}

std::string throwMsgOnRowsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    std::ostringstream os;
    os << matName << " should have " << std::to_string(ps->uDim) << " or " << std::to_string(ps->fullUDim)
       << " rows. But " << matName << " has " << std::to_string(mat.rows()) << " rows.";
    return os.str();
}

std::string throwMsgOnColsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    std::ostringstream os;
    os << matName << " should have " << std::to_string(ps->xDim) << " or " << std::to_string(ps->fullXDim)
       << " cols. But " << matName << " has " << std::to_string(mat.cols()) << " cols.";
    return os.str();
}

std::string throwMsgOnColsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    std::ostringstream os;
    os << matName << " should have " << std::to_string(ps->uDim) << " or " << std::to_string(ps->fullUDim)
       << " cols. But " << matName << " has " << std::to_string(mat.cols()) << " cols.";
    return os.str();
}

std::string throwMsgOnColsOnPSxuDim(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2, const PreviewSystem* ps)
{
    std::ostringstream os;
    os << mat1Name << " and " << mat2Name << " should respectively have "
       << std::to_string(ps->xDim) << " and " << std::to_string(ps->uDim) << " cols. But "
       << mat1Name << " has " << std::to_string(mat1.cols()) << " cols and "
       << mat2Name << " has " << std::to_string(mat2.cols()) << " cols.";
    return os.str();
}

std::string throwMsgOnColsOnPSXUDim(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2, const PreviewSystem* ps)
{
    std::ostringstream os;
    os << mat1Name << " and " << mat2Name << " should respectively have "
       << std::to_string(ps->xDim) << " and " << std::to_string(ps->uDim) << " cols or "
       << std::to_string(ps->fullXDim) << " and " << std::to_string(ps->fullUDim) << " cols. But "
       << mat1Name << " has " << std::to_string(mat1.cols()) << " cols and "
       << mat2Name << " has " << std::to_string(mat2.cols()) << " cols.";
    return os.str();
}

std::string throwMsgOnBadNewDim(const Eigen::MatrixXd& mat, Eigen::Index new_dim)
{
    std::ostringstream os;
    os << "Can not span the matrix to the desired new dimension. The  desired dimension ("
       << new_dim << ") is not a multiple of the number of rows ("
       << mat.rows() << ").";
    return os.str();
}

} // namespace copra
