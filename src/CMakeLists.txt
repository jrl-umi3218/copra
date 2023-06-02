#
# Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

set(SOURCES
    AutoSpan.cpp
    costFunctions.cpp
    PreviewSystem.cpp
    QuadProgSolver.cpp
    SolverInterface.cpp
    constraints.cpp
    debugUtils.cpp
    LMPC.cpp
    InitialStateLMPC.cpp
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
    ../include/InitialStateLMPC.h
    ../include/solverUtils.h
    ../include/typedefs.h)

set(COPRA_COMPILE_DEFS "")
if(${eigen-lssol_FOUND})
  list(APPEND COPRA_COMPILE_DEFS EIGEN_LSSOL_FOUND)
  list(APPEND SOURCES LSSOLSolver.cpp)
  list(APPEND HEADERS ../include/LSSOLSolver.h)
endif()

if(${eigen-gurobi_FOUND})
  list(APPEND COPRA_COMPILE_DEFS EIGEN_GUROBI_FOUND)
  list(APPEND SOURCES GUROBISolver.cpp)
  list(APPEND HEADERS ../include/GUROBISolver.h)
endif()

if(${eigen-qld_FOUND})
  list(APPEND COPRA_COMPILE_DEFS EIGEN_QLD_FOUND)
  list(APPEND SOURCES QLDSolver.cpp)
  list(APPEND HEADERS ../include/QLDSolver.h)
endif()

if(${eigen-osqp_FOUND})
  list(APPEND COPRA_COMPILE_DEFS EIGEN_OSQP_FOUND)
  list(APPEND SOURCES OSQPSolver.cpp)
  list(APPEND HEADERS ../include/OSQPSolver.h)
endif()

add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})
if(MSVC)
  target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX)
endif()
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS
                                                 "-Dcopra_EXPORTS")
set_target_properties(
  ${PROJECT_NAME} PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR}
                             VERSION ${PROJECT_VERSION})
target_compile_definitions(${PROJECT_NAME}
                           PRIVATE $<BUILD_INTERFACE:EIGEN_MPL2_ONLY>)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${COPRA_COMPILE_DEFS})
target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} PUBLIC eigen-quadprog::eigen-quadprog)
link_optional_qps(${PROJECT_NAME})
if(TARGET Eigen3::Eigen)
  target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
else()
  target_include_directories(${PROJECT_NAME} SYSTEM
                             PUBLIC "${EIGEN3_INCLUDE_DIR}")
endif()

install(
  TARGETS ${PROJECT_NAME}
  EXPORT "${TARGETS_EXPORT_NAME}"
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DESTINATION})
