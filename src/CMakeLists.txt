#To complete
set(SOURCES SolverInterface.cpp QuadProgSolver.cpp QLDSolver.cpp previewController.cpp)
set(HEADERS SolverInterface.h QuadProgSolver.h QLDSolver.h previewController.h solverUtils.h)

if(${EIGEN_LSSOL_FOUND})
    list(APPEND SOURCES LSSOLSolver.cpp)
    list(APPEND HEADERS LSSOLSolver.h)
endif()

include_directories(${Boost_INCLUDE_DIRS})
add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})

target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-qld)
PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-quadprog)
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-lssol)
    add_definitions(-DLSSOL_SOLVER_FOUND)
endif()

set(INSTALL_PATH include/PC)
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS ${PROJECT_NAME} DESTINATION "lib")