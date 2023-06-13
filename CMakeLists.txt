#
# Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

# Version minimum
cmake_minimum_required(VERSION 3.1.3)

set(PROJECT_NAME copra)
set(PROJECT_DESCRIPTION
    "This library is an implementation of a linear model predictive control")
set(PROJECT_URL "https://github.com/jrl-umi3218/copra")
set(PROJECT_DEBUG_POSTFIX "_d")
set(PROJECT_VERSION 1.3.1)

include(cmake/base.cmake)
include(cmake/msvc-specific.cmake)
include(optional_qps.cmake)

# SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")
set(CMAKE_CXX_STANDARD 14)

option(PYTHON_BINDING "Generate python bindings." OFF)
option(BUILD_TESTING "Build unit tests." ON)

project(
  ${PROJECT_NAME}
  LANGUAGES CXX
  VERSION ${PROJECT_VERSION})

set_default_cmake_build_type(RelWithDebInfo)

set(INSTALL_GENERATED_HEADERS OFF)

if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

set(BOOST_REQUIRED 1.48)
set(Eigen_REQUIRED "eigen3 >= 3.2.0")

add_project_dependency(eigen-quadprog REQUIRED)
add_optional_qps(COPRA_PACKAGE_LIST)
foreach(COPRA_PKG ${COPRA_PACKAGE_LIST})
  list(APPEND _PACKAGE_CONFIG_DEPENDENCIES
       "find_package(${COPRA_PKG} REQUIRED)")
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
