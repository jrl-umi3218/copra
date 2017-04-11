# Version minimum
cmake_minimum_required(VERSION 2.8)
include(cmake/base.cmake)
include(cmake/eigen.cmake)
include(cmake/boost.cmake)

set(PROJECT_NAME mpc)
set(PROJECT_DESCRIPTION "This a generic model preview control")
set(PROJECT_URL "https://github.com/vsamy/preview_controller")

#SET(CXX_DISABLE_WERROR True)
set(DOXYGEN_USE_MATHJAX "YES")

# Check compiler version
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0")
        message(FATAL_ERROR "Insufficient gcc version (version superior to 5.0 is required)")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "3.4")
        message(FATAL_ERROR "Insufficient clang version (version superior to 3.4 is required)")
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "14.0")
        message(FATAL_ERROR "Insufficient msvc version (version superior to 14.0 is required)")
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

if(${_Eigen_VERSION} VERSION_LESS 3.3)
  add_definitions(-DEIGEN_3_2_FOUND) # Will add Eigen::Index as in 3.3
endif()

ADD_REQUIRED_DEPENDENCY("eigen-quadprog")
ADD_OPTIONAL_DEPENDENCY("eigen-qld")
ADD_OPTIONAL_DEPENDENCY("eigen-lssol")
ADD_OPTIONAL_DEPENDENCY("eigen-gurobi")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++14 -pedantic")

# Header for define and include found c++ QP solver
configure_file(cmake/solverConfig.cmake.in include/${HEADER_DIR}/solverConfig.h)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/${HEADER_DIR}/solverConfig.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${HEADER_DIR})

# Include the automatically generated configuration files.
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include/${HEADER_DIR})

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(binding/python)

# Add dependency towards the library in the pkg-config file.
pkg_config_append_libs(${PROJECT_NAME})

SETUP_PROJECT_FINALIZE()