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

# Version minimum
cmake_minimum_required(VERSION 3.12.4) # Facilitate everything for python and boost 
include(cmake/base.cmake)
include(cmake/msvc-specific.cmake)
include(optional_qps.cmake)

set(PROJECT_NAME copra)
set(PROJECT_DESCRIPTION "This library is an implementation of a linear model predictive control")
set(PROJECT_URL "https://github.com/vsamy/copra")
set(PROJECT_DEBUG_POSTFIX "_d")

# SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")
set(CMAKE_CXX_STANDARD 14)

option(PYTHON_BINDING "Generate python bindings." OFF)
option(BUILD_TESTING "Build unit tests." ON)

project(${PROJECT_NAME} CXX)

setup_project()
set(INSTALL_GENERATED_HEADERS OFF)

if(MSVC)
    set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
    # add_compile_options("/bigobj")
endif(MSVC)

set(BOOST_REQUIRED 1.48)
set(Eigen_REQUIRED "eigen3 >= 3.2.0")

set(COPRA_PACKAGE_LIST Eigen3 eigen-quadprog)
add_optional_qps(${COPRA_PACKAGE_LIST})
foreach(COPRA_PKG ${COPRA_PACKAGE_LIST})
    add_project_dependency(${COPRA_PKG} REQUIRED)
endforeach()

add_subdirectory(src)

if(${BUILD_TESTING})
    add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
    add_subdirectory(binding/python)
endif()

setup_project_finalize()
setup_project_package_finalize()
