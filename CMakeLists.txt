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
cmake_minimum_required(VERSION 2.8)
include(cmake/base.cmake)
include(cmake/eigen.cmake)
include(cmake/boost.cmake)
include(cmake/msvc-specific.cmake)

set(PROJECT_NAME copra)
set(PROJECT_DESCRIPTION "This library is an implementation of a linear model predictive control")
set(PROJECT_URL "https://github.com/vsamy/copra")
set(PROJECT_DEBUG_POSTFIX "_d")

# SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

cmake_policy(VERSION ${CMAKE_VERSION})

project(${PROJECT_NAME} CXX)
set (CMAKE_CXX_STANDARD 11)


# # Check compiler version
# if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#     if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "5.0")
#         message(FATAL_ERROR "Insufficient gcc version (version superior to 5.0 is required)")
#     endif()
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
#     if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "3.4")
#         message(FATAL_ERROR "Insufficient clang version (version superior to 3.4 is required)")
#     endif()
# elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
#     if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "14.0")
#         message(FATAL_ERROR "Insufficient msvc version (version superior to 14.0 is required)")
#     endif()
# endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" lower_case)
if(lower_case STREQUAL "debug")
    set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG")
endif()

option(PYTHON_BINDING "Generate python bindings." OFF)
option(DISABLE_TESTS "Build unit tests." OFF)

setup_project()

if(MSVC)
    set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
    add_compile_options("/bigobj")
endif(MSVC)

set(CMAKE_CXX_STANDARD 11)

set(BOOST_REQUIRED 1.48)
set(BOOST_COMPONENTS system timer)
set(Eigen_REQUIRED "eigen3 >= 3.2.0")

#########################
# External dependencies
if(PKG_CONFIG_FOUND)
  search_for_eigen()
else()
  find_path(Eigen_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
    PATHS
    ${CMAKE_INSTALL_PREFIX}/include
    PATH_SUFFIXES eigen3 eigen
  )
  include_directories(SYSTEM ${Eigen_INCLUDE_DIR})
endif()
#########################

add_required_dependency("eigen-quadprog" "eigen-quadprog_d")
add_optional_dependency("eigen-qld" "eigen-qld_d")
add_optional_dependency("eigen-lssol" "eigen-lssol_d")
add_optional_dependency("eigen-gurobi" "eigen-gurobi_d")
if(${EIGEN_LSSOL_FOUND})
    add_compile_options(-DEIGEN_LSSOL_FOUND)
    pkg_config_append_cflags(-DEIGEN_LSSOL_FOUND)
endif()
if(${EIGEN_GUROBI_FOUND})
    add_compile_options(-DEIGEN_GUROBI_FOUND)
    pkg_config_append_cflags(-DEIGEN_GUROBI_FOUND)
    endif()
if(${EIGEN_QLD_FOUND})
    add_compile_options(-DEIGEN_QLD_FOUND)
    pkg_config_append_cflags(-DEIGEN_QLD_FOUND)
endif()

# if(NOT WIN32)
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -pedantic")
# endif()

add_subdirectory(src)
if(NOT ${DISABLE_TESTS})
    add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
    add_subdirectory(binding/python)
endif()

# Add dependency towards the library in the pkg-config file.
pkg_config_append_libs(${PROJECT_NAME})

setup_project_finalize()
