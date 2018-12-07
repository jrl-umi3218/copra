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

#include <Eigen/Core>

namespace {

Eigen::MatrixXd spanMatrix(const Eigen::MatrixXd& m, int size, int addCols = 0)
{
    Eigen::MatrixXd mout = Eigen::MatrixXd::Zero(m.rows() * size, m.cols() * (size + addCols));
    for (int i = 0; i < size; ++i)
        mout.block(i * m.rows(), i * m.cols(), m.rows(), m.cols()) = m;
    return mout;
}

Eigen::VectorXd spanVector(const Eigen::VectorXd& v, int size)
{
    Eigen::VectorXd vout = Eigen::VectorXd::Zero(v.rows() * size);
    for (int i = 0; i < size; ++i)
        vout.segment(i * v.rows(), v.rows()) = v;
    return vout;
}

} // namespace