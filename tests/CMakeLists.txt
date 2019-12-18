#
# Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

macro(addTest testName)
    add_executable(${testName} ${testName}.cpp systems.h tools.h tools.cpp)
    target_link_libraries(${testName} PRIVATE Boost::unit_test_framework ${PROJECT_NAME})
    add_test(${testName}Unit ${testName})
    # Adding a project configuration file (for MSVC only)
    generate_msvc_dot_user_file(${testName} "$(SolutionDir)/src/${PROJECT_NAME}")
endmacro(addTest)

addTest(TestSolvers)
addTest(TestLMPC)
