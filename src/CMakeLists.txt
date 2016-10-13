#To complete
set(SOURCES SolverInterface.cpp QuadProgSolver.cpp QLDSolver.cpp)
set(HEADERS SolverInterface.h QuadProgSolver.h QLDSolver.h)

if(${EIGEN_LSSOL_FOUND})
    list(APPEND SOURCES LSSOLSolver.cpp)
    list(APPEND HEADERS LSSOLSolver.h)
endif()

add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})

PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-qld)
PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-quadprog)
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-lssol)
    add_definitions(-DLSSOL_SOLVER_FOUND)
endif()

set(INSTALL_PATH include/PC)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS ${PROJECT_NAME} DESTINATION "lib")