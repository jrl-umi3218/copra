# Version minimum
cmake_minimum_required(VERSION 2.8)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/boost.cmake)

SET(PROJECT_NAME PreviewControl)
SET(PROJECT_DESCRIPTION "...")
SET(PROJECT_URL "")

#SET(CXX_DISABLE_WERROR True)
SET(DOXYGEN_USE_MATHJAX "YES")

project(${PROJECT_NAME} CXX)
SETUP_PROJECT()

SEARCH_FOR_EIGEN()
ADD_REQUIRED_DEPENDENCY("eigen-quadprog")
ADD_REQUIRED_DEPENDENCY("eigen-qld")
ADD_OPTIONAL_DEPENDENCY("eigen-lssol")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1y -pedantic")

add_subdirectory(src)
add_subdirectory(tests)

# untill the lib works fine
SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/_build/tmp_install)
SET(CMAKE_BUILD_TYPE DEBUG) 

SETUP_PROJECT_FINALIZE()