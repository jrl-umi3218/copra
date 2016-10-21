#To complete
set(SOURCES SolverInterface.cpp QuadProgSolver.cpp QLDSolver.cpp 
    PreviewController.cpp PreviewSystem.cpp Constrains.cpp)
set(HEADERS SolverInterface.h QuadProgSolver.h QLDSolver.h 
    PreviewController.h PreviewSystem.h Constrains.h
    solverUtils.h)

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

set(INSTALL_PATH include/${PROJECT_NAME})
install(FILES ${HEADERS} DESTINATION ${INSTALL_PATH})
install(TARGETS ${PROJECT_NAME} DESTINATION "lib")