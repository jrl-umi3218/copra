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

set(SOURCES
    AutoSpan.cpp
    costFunctions.cpp
    PreviewSystem.cpp
    QuadProgSolver.cpp
    SolverInterface.cpp
    constraints.cpp
    debugUtils.cpp
    previewController.cpp
    solverUtils.cpp)

set(HEADERS 
    ../include/AutoSpan.h
    ../include/costFunctions.h
    ../include/PreviewSystem.h 
    ../include/QuadProgSolver.h 
    ../include/SolverInterface.h 
    ../include/constraints.h
    ../include/debugUtils.h
    ../include/previewController.h 
    ../include/solverUtils.h 
    ../include/typedefs.h)

if(${EIGEN_LSSOL_FOUND})
    list(APPEND SOURCES LSSOLSolver.cpp)
    list(APPEND HEADERS ../include/LSSOLSolver.h)
endif()

if(${EIGEN_GUROBI_FOUND})
    list(APPEND SOURCES GUROBISolver.cpp)
    list(APPEND HEADERS ../include/GUROBISolver.h)
endif()

if(${EIGEN_QLD_FOUND})
    list(APPEND SOURCES QLDSolver.cpp)
    list(APPEND HEADERS ../include/QLDSolver.h)
endif()

add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})

PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-quadprog)
if(${EIGEN_QLD_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-qld)
endif()
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-lssol)
endif()
if(${EIGEN_GUROBI_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-gurobi)
endif()

set(INSTALL_PATH include/${PROJECT_NAME})
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)