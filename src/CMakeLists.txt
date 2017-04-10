#To complete
set(SOURCES
    PreviewSystem.cpp
    QuadProgSolver.cpp
    SolverInterface.cpp
    constraints.cpp
    debugUtils.cpp
    previewController.cpp
    solverUtils.cpp)

set(HEADERS 
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

include_directories(${Boost_INCLUDE_DIRS})
add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
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