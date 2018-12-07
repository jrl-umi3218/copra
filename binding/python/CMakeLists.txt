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

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/pyCopra/__init__.py COPYONLY)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-std=c++14") # g++ fails with -Werror on boost_python for whatever reason...
endif()

findpython()
set(BOOST_COMPONENTS system numpy)
search_for_boost()

if(NOT PYTHON_INSTALL_PATH)
    set(PYTHON_INSTALL_PATH ${PYTHON_SITELIB})
endif()

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(pyCopra SHARED CopraBindings.cpp)
target_link_libraries(pyCopra ${PROJECT_NAME} ${PYTHON_LIBRARY} ${Boost_PYTHON_LIBRARIES})
pkg_config_use_dependency(pyCopra eigen-quadprog)
if(${EIGEN_QLD_FOUND})
    pkg_config_use_dependency(pyCopra eigen-qld)
    add_compile_options(-DEIGEN_QLD_FOUND)
endif()
if(${EIGEN_LSSOL_FOUND})
    pkg_config_use_dependency(pyCopra eigen-lssol)
    add_compile_options(-DEIGEN_LSSOL_FOUND)
endif()
if(${EIGEN_GUROBI_FOUND})
    pkg_config_use_dependency(pyCopra eigen-gurobi)
    add_compile_options(-DEIGEN_GUROBI_FOUND)
endif()
if(WIN32)
    set_target_properties(pyCopra PROPERTIES SUFFIX ".pyd")
else()
    set_target_properties(pyCopra PROPERTIES PREFIX "")
endif()

install(TARGETS pyCopra DESTINATION "${PYTHON_INSTALL_PATH}/pyCopra")
python_install_build(pyCopra __init__.py "${PYTHON_INSTALL_PATH}")

if(NOT ${DISABLE_TESTS})
    add_subdirectory(tests)
endif()