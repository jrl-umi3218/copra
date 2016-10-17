enable_testing()

find_package(Boost COMPONENTS unit_test_framework REQUIRED)

add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

add_executable(TestSolvers TestSolvers.cpp )
target_link_libraries(TestSolvers ${Boost_LIBRARIES} ${PROJECT_NAME})
PKG_CONFIG_USE_DEPENDENCY(TestSolvers eigen-qld)
PKG_CONFIG_USE_DEPENDENCY(TestSolvers eigen-quadprog)
if(${EIGEN_LSSOL_FOUND})
    PKG_CONFIG_USE_DEPENDENCY(TestSolvers eigen-lssol)
    add_definitions(-DLSSOL_SOLVER_FOUND)
endif()

add_test(TestSolversUnit TestSolvers)