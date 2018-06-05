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

namespace Eigen {

#if EIGEN_VERSION_AT_LEAST(3, 2, 0)
using Index = Eigen::Matrix<int, 1, 1>::Index;
#endif

} // namespace Eigen

namespace copra {

template <typename T1, typename T2 = std::true_type, typename T3 = std::true_type>
struct is_all_arithmetic : is_all_arithmetic<typename std::is_arithmetic<std::decay_t<T1> >::type, T2, T3> {
};

template <typename T2, typename T3>
struct is_all_arithmetic<std::false_type, T2, T3> {
    static const bool value = false;
};

template <typename T2, typename T3>
struct is_all_arithmetic<std::true_type, T2, T3> : is_all_arithmetic<typename std::is_arithmetic<std::decay_t<T2> >::type, T3, std::true_type> {
};

template <>
struct is_all_arithmetic<std::true_type, std::true_type, std::true_type> {
    static const bool value = true;
};

} // namespace copra
