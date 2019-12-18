/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#include "AutoSpan.h"
#include "debugUtils.h"

namespace copra {

void AutoSpan::spanMatrix(Eigen::MatrixXd& mat, Eigen::Index new_dim, int addCols)
{
    auto matRows = mat.rows();
    if (new_dim == matRows)
        return;

    auto matCols = mat.cols();
    auto tmp = mat;
    auto nrStep = new_dim / matRows;
    if (nrStep * matRows != new_dim)
        DOMAIN_ERROR_EXCEPTION(throwMsgOnBadNewDim(mat, new_dim));

    mat = Eigen::MatrixXd::Zero(new_dim, matCols * (nrStep + addCols));
    for (auto i = 0; i < nrStep; ++i)
        mat.block(i * matRows, i * matCols, matRows, matCols) = tmp;
}

void AutoSpan::spanVector(Eigen::VectorXd& vec, Eigen::Index new_dim)
{
    auto vecRows = vec.rows();
    if (new_dim == vecRows)
        return;

    auto nrStep = new_dim / vecRows;
    auto tmp = vec;
    if (nrStep * vecRows != new_dim)
        DOMAIN_ERROR_EXCEPTION(throwMsgOnBadNewDim(vec, new_dim));

    vec.conservativeResize(new_dim);
    for (auto i = 1; i < nrStep; ++i)
        vec.segment(i * vecRows, vecRows) = tmp;
}

} // namespace copra
