include(../../cmake/python.cmake)

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/mpcontroller/__init__.py COPYONLY)

find_package(Boost COMPONENTS python timer REQUIRED)
FINDPYTHON()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(_mpcontroller SHARED MPController.cpp)
target_link_libraries(_mpcontroller ${Boost_LIBRARIES} ${PROJECT_NAME})
PKG_CONFIG_USE_DEPENDENCY(_mpcontroller eigen-qld)
PKG_CONFIG_USE_DEPENDENCY(_mpcontroller eigen-quadprog)
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(_mpcontroller eigen-lssol)
    add_definitions(-DLSSOL_SOLVER_FOUND)
endif()
set_target_properties(_mpcontroller PROPERTIES PREFIX "")

install(TARGETS _mpcontroller DESTINATION "${PYTHON_SITELIB}/mpcontroller")
PYTHON_INSTALL_BUILD(mpcontroller __init__.py "${PYTHON_SITELIB}")