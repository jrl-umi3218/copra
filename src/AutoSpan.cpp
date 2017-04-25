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

// header
#include "AutoSpan.h"

namespace mpc {

void AutoSpan::spanMatrix(Eigen::MatrixXd &mat, Eigen::Index new_dim, int addCols) const
{
    auto matRows = mat.rows();
    if (new_dim == matRows)
        return;

    auto matCols = mat.cols();
    auto tmp = mat;
    auto nrStep = max_dim / matRows;
    if (nrStep * matRows != max_dim)
        DOMAIN_ERROR_EXCEPTION(mat, new_dim);

    mat = Eigen::MatrixXd::Zero(new_dim, matCols * (nrStep + addCols));
    for (auto i = 0; i < nrStep; ++i)
        mat.block(i * matRows, i * matCols, matRows, matCols) = tmp;
}

void AutoSpan::spanVector(Eigen::VectorXd &vec, Eigen::Index new_dim) const
{
    auto vecRows = vec.rows();
    if (max_dim == vecRows)
        return;

    auto nrStep = new_dim / vecRows;
    auto tmp = vec;
    if (nrStep * vecRows != max_dim)
        DOMAIN_ERROR_EXCEPTION(vec, new_dim);

    vec.conservativeResize(new_dim);
    for (auto i = 1; i < nrStep; ++i)
        vec.segment(i * vecRows, vecRows) = tmp;
}

} // namespace mpc