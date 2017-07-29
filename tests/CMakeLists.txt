# This file is part of copra.

# copra is free software: you can redistribute it and/or
# modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# copra is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with copra.  If not, see
# <http://www.gnu.org/licenses/>.

enable_testing()

set(BOOST_COMPONENTS system unit_test_framework)
SEARCH_FOR_BOOST()

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

macro(addTest testName)
    add_executable(${testName} ${testName}.cpp systems.h tools.h)
    target_link_libraries(${testName} ${Boost_LIBRARIES} ${PROJECT_NAME})
    PKG_CONFIG_USE_DEPENDENCY(${testName} eigen-quadprog)
    if(${EIGEN_QLD_FOUND})
        PKG_CONFIG_USE_DEPENDENCY(${testName} eigen-qld)
    endif()
    if(${EIGEN_LSSOL_FOUND})
        PKG_CONFIG_USE_DEPENDENCY(${testName} eigen-lssol)
    endif()
    if(${EIGEN_GUROBI_FOUND})
        PKG_CONFIG_USE_DEPENDENCY(${testName} eigen-gurobi)
    endif()
    add_test(${testName}Unit ${testName})
endmacro(addTest)

addTest(TestSolvers)
addTest(TestLMPC)