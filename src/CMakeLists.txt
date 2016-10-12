#To complete
set(SOURCES solvers.cpp)
set(HEADERS solvers.h)

add_library(${PROJECT_NAME} SHARED ${SOURCES} ${HEADERS})
PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-qld)
PKG_CONFIG_USE_DEPENDENCY(${PROJECT_NAME} eigen-quadprog)