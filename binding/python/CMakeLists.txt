# This file is part of mpc.

# mpc is free software: you can redistribute it and/or
# modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mpc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with mpc.  If not, see
# <http://www.gnu.org/licenses/>.

include(../../cmake/python.cmake)

SEARCH_FOR_BOOST()

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/${HEADER_DIR}/__init__.py COPYONLY)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-std=c++14") # g++ fails with -Werror on boost_python for whatever reason...
endif()

find_package(Boost COMPONENTS python REQUIRED)
FINDPYTHON()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(_mpc SHARED MPController.cpp)
target_link_libraries(_mpc ${Boost_LIBRARIES} ${PROJECT_NAME})
PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-quadprog)
if(${EIGEN_QLD_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-qld)
endif()
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-lssol)
endif()
if(${EIGEN_GUROBI_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_mpc eigen-gurobi)
endif()
set_target_properties(_mpc PROPERTIES PREFIX "")

install(TARGETS _mpc DESTINATION "${PYTHON_SITELIB}/${HEADER_DIR}")
PYTHON_INSTALL_BUILD(mpc __init__.py "${PYTHON_SITELIB}")