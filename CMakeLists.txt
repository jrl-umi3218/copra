# Version minimum
cmake_minimum_required(VERSION 2.8)
include(cmake/base.cmake)
include(cmake/eigen.cmake)
include(cmake/boost.cmake)

set(PROJECT_NAME PreviewControl)
set(PROJECT_DESCRIPTION "...")
set(PROJECT_URL "")

#SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

project(${PROJECT_NAME} CXX)
SETUP_PROJECT()

set(BOOST_REQUIRED 1.48)
set(BOOST_COMPONENTS system timer)
SEARCH_FOR_BOOST()
SEARCH_FOR_EIGEN()

ADD_REQUIRED_DEPENDENCY("eigen-quadprog")
ADD_REQUIRED_DEPENDENCY("eigen-qld")
ADD_OPTIONAL_DEPENDENCY("eigen-lssol")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1y -pedantic")

add_subdirectory(src)
add_subdirectory(tests)

# untill the lib works fine
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/_build/tmp_install)
set(CMAKE_BUILD_TYPE DEBUG) 

SETUP_PROJECT_FINALIZE()