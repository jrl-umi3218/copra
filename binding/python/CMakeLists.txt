include(../../cmake/python.cmake)

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/mpc/__init__.py COPYONLY)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-std=c++14") # g++ fails with -Werror on boost_python for whatever reason...
endif()

find_package(Boost COMPONENTS python timer REQUIRED)
FINDPYTHON()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(_mpc SHARED MPController.cpp)
target_link_libraries(_mpc ${Boost_LIBRARIES} ${PROJECT_NAME})
PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-quadprog)
if(${EIGEN_QLD_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-qld)
    add_definitions(-DQLD_SOLVER_FOUND)
endif()
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-lssol)
    add_definitions(-DLSSOL_SOLVER_FOUND)
endif()
if(${GUROBI_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-gurobi)
    add_definitions(-DGUROBI_SOLVER_FOUND)
endif()
set_target_properties(_mpc PROPERTIES PREFIX "")

install(TARGETS _mpc DESTINATION "${PYTHON_SITELIB}/mpc")
PYTHON_INSTALL_BUILD(mpc __init__.py "${PYTHON_SITELIB}")