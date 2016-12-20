# Version minimum
cmake_minimum_required(VERSION 2.8)
include(cmake/base.cmake)
include(cmake/eigen.cmake)
include(cmake/boost.cmake)

set(PROJECT_NAME modelpreviewcontrol)
set(PROJECT_DESCRIPTION "This a generic model preview control")
set(PROJECT_URL "https://github.com/vsamy/preview_controller")

#SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

# Check compiler version
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0")
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.4")
        message(FATAL_ERROR "Insufficient clang version")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "14.0")
        message(FATAL_ERROR "Insufficient msvc version")
    endif()
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" lower_case)
if(lower_case STREQUAL "debug")
    set(CMAKE_CXX_FLAGS_DEBUG "-D_DEBUG")
endif()

project(${PROJECT_NAME} CXX)
SETUP_PROJECT()

set(BOOST_REQUIRED 1.48)
set(BOOST_COMPONENTS system timer)
set(Eigen_REQUIRED "eigen3 >= 3.2.7") # This eigen version allows move semantics.
SEARCH_FOR_BOOST()
SEARCH_FOR_EIGEN()

ADD_REQUIRED_DEPENDENCY("eigen-quadprog")
ADD_OPTIONAL_DEPENDENCY("eigen-qld")
ADD_OPTIONAL_DEPENDENCY("eigen-lssol")
ADD_OPTIONAL_DEPENDENCY("eigen-gurobi")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -pedantic")

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(binding/python)

SETUP_PROJECT_FINALIZE()