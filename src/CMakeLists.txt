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

set(SOURCES
    AutoSpan.cpp
    costFunctions.cpp
    PreviewSystem.cpp
    QuadProgSolver.cpp
    SolverInterface.cpp
    constraints.cpp
    debugUtils.cpp
    LMPC.cpp
    solverUtils.cpp)

set(HEADERS 
    ../include/api.h
    ../include/AutoSpan.h
    ../include/costFunctions.h
    ../include/PreviewSystem.h 
    ../include/QuadProgSolver.h 
    ../include/SolverInterface.h 
    ../include/constraints.h
    ../include/debugUtils.h
    ../include/LMPC.h 
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
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-Dcopra_EXPORTS")
target_include_directories(${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PUBLIC eigen-quadprog)
link_optional_qps(${PROJECT_NAME})
if(TARGET Eigen3::Eigen) 
  target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen) 
else() 
  target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}") 
endif() 

set(INSTALL_PATH include/${PROJECT_NAME})
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)