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
#include "debugUtils.h"

// stl
#include <sstream>

namespace mpc {

void checkMat(const char* isMatName, const Eigen::MatrixXd& isMat, const Eigen::MatrixXd& shouldBeMat)
{
    if (isMat.rows() != shouldBeMat.rows() || isMat.cols() != shouldBeMat.cols()) {
        std::ostringstream os;
        os << "Bad dimension for " << isMatName
           << ". It should be an (" << std::to_string(shouldBeMat.rows()) << "-by-" << std::to_string(shouldBeMat.cols())
           << ") matrix but you gave an (" << std::to_string(isMat.rows()) << "-by-" << std::to_string(isMat.cols()) << ") matrix";
        throw std::domain_error(os.str());
    }
}

void checkRows(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2)
{
    if (mat1.rows() != mat2.rows()) {
        std::ostringstream os;
        os << mat1Name << " and " << mat2Name << " should have same number of rows. "
           << mat1Name << " has " << std::to_string(mat1.rows()) << " rows and "
           << mat2Name << " has " << std::to_string(mat2.rows()) << " rows.";
        throw std::domain_error(os.str());
    }
}

void checkSquaredMat(const char* matName, const Eigen::MatrixXd& mat)
{
    if (mat.rows() != mat.cols()) {
        std::ostringstream os;
        os << matName << " should be a squared matrix or it is a ("
           << mat.rows() << "-by-"
           << mat.cols() << ") matrix.";
        throw std::domain_error(os.str());
    }
}

void checkRowsOnDim(const char* matName, const Eigen::MatrixXd& mat, Eigen::Index dim)
{
    if (mat.rows() != dim) {
        std::ostringstream os;
        os << matName << " should have exactly " << std::to_string(dim) << " rows. But "
           << matName << " has " << std::to_string(mat.rows()) << " rows.";
        throw std::domain_error(os.str());
    }
}

void checkRowsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    if (mat.rows() != ps->xDim && mat.rows() != ps->fullXDim) {
        std::ostringstream os;
        os << matName << " should have " << std::to_string(ps->xDim) << " or " << std::to_string(ps->fullXDim)
           << " rows. But it has " << matName << " has " << std::to_string(mat.rows()) << " rows.";
        throw std::domain_error(os.str());
    }
}

void checkRowsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    if (mat.rows() != ps->uDim && mat.rows() != ps->fullUDim) {
        std::ostringstream os;
        os << matName << " should have " << std::to_string(ps->uDim) << " or " << std::to_string(ps->fullUDim)
           << " rows. But it has " << matName << " has " << std::to_string(mat.rows()) << " rows.";
        throw std::domain_error(os.str());
    }
}

void checkColsOnPSXDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    if (mat.cols() != ps->xDim && mat.cols() != ps->fullXDim) {
        std::ostringstream os;
        os << matName << " should have " << std::to_string(ps->xDim) << " or " << std::to_string(ps->fullXDim)
           << " cols. But it has " << matName << " has " << std::to_string(mat.cols()) << " cols.";
        throw std::domain_error(os.str());
    }
}

void checkColsOnPSUDim(const char* matName, const Eigen::MatrixXd& mat, const PreviewSystem* ps)
{
    if (mat.cols() != ps->uDim && mat.cols() != ps->fullUDim) {
        std::ostringstream os;
        os << matName << " should have " << std::to_string(ps->uDim) << " or " << std::to_string(ps->fullUDim)
           << " cols. But it has " << matName << " has " << std::to_string(mat.cols()) << " cols.";
        throw std::domain_error(os.str());
    }
}

void checkColsOnPSXUDim(const char* mat1Name, const char* mat2Name, const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2, const PreviewSystem* ps)
{
    if ((mat1.cols() != ps->xDim && mat2.cols() != ps->uDim) || (mat1.cols() != ps->fullXDim && mat2.cols() != ps->fullUDim)) {
        std::ostringstream os;
        os << mat1Name << " and " << mat2Name << " should respectively have "
           << std::to_string(ps->xDim) << " and " << std::to_string(ps->uDim) << " cols or "
           << std::to_string(ps->fullXDim) << " and " << std::to_string(ps->fullUDim) << " cols. But "
           << mat1Name << " has " << std::to_string(mat1.cols()) << " cols and "
           << mat2Name << " has " << std::to_string(mat2.cols()) << "cols.";
        throw std::domain_error(os.str());
    }
}
}