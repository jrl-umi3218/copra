#
# Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

macro(addTest testName)
    add_executable(${testName} ${testName}.cpp systems.h tools.h tools.cpp doctest.h)
    target_link_libraries(${testName} PUBLIC ${PROJECT_NAME})
    target_compile_definitions(${testName} PRIVATE DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN)
    add_test(${testName}Unit ${testName})
    # Adding a project configuration file (for MSVC only)
    generate_msvc_dot_user_file(${testName} "$(SolutionDir)/src/${PROJECT_NAME}")
endmacro(addTest)

addTest(TestSolvers)
addTest(TestLMPC)
