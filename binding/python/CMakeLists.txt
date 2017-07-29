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

include(../../cmake/python.cmake)

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/${HEADER_DIR}/__init__.py COPYONLY)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-std=c++14") # g++ fails with -Werror on boost_python for whatever reason...
endif()

set(BOOST_COMPONENTS system python)
SEARCH_FOR_BOOST()
FINDPYTHON()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(_copra SHARED CopraBindings.cpp)
target_link_libraries(_copra ${Boost_LIBRARIES} ${PROJECT_NAME})
PKG_CONFIG_USE_DEPENDENCY(_copra eigen-quadprog)
if(${EIGEN_QLD_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_copra eigen-qld)
endif()
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_copra eigen-lssol)
endif()
if(${EIGEN_GUROBI_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_copra eigen-gurobi)
endif()
set_target_properties(_copra PROPERTIES PREFIX "")

install(TARGETS _copra DESTINATION "${PYTHON_SITELIB}/${HEADER_DIR}")
PYTHON_INSTALL_BUILD(copra __init__.py "${PYTHON_SITELIB}")