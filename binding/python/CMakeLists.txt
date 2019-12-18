#
# Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/pyCopra/__init__.py COPYONLY)

set(Boost_USE_STATIC_LIBS ON)
find_package(Python REQUIRED COMPONENTS Development)
find_package(Boost REQUIRED COMPONENTS python${Python_VERSION_MAJOR}${Python_VERSION_MINOR} numpy)

if(NOT PYTHON_INSTALL_PATH)
    set(PYTHON_INSTALL_PATH ${PYTHON_SITELIB})
endif()

add_library(pyCopra SHARED CopraBindings.cpp)
target_compile_definitions(pyCopra PUBLIC -DBOOST_NUMPY_STATIC_LIB -DBOOST_PYTHON_STATIC_LIB)
target_link_libraries(pyCopra PUBLIC eigen-quadprog::eigen-quadprog Python::Python Boost::python${Python_VERSION_MAJOR}${Python_VERSION_MINOR} Boost::numpy ${PROJECT_NAME})
link_optional_qps(pyCopra)
if(TARGET Eigen3::Eigen) 
  target_link_libraries(pyCopra PUBLIC Eigen3::Eigen) 
else() 
  target_include_directories(pyCopra SYSTEM PUBLIC "${EIGEN3_INCLUDE_DIR}") 
endif() 
if(WIN32)
    set_target_properties(pyCopra PROPERTIES SUFFIX ".pyd")
else()
    set_target_properties(pyCopra PROPERTIES PREFIX "")
endif()

install(TARGETS pyCopra DESTINATION "${PYTHON_INSTALL_PATH}/pyCopra")
install(FILES __init__.py DESTINATION "${PYTHON_INSTALL_PATH}/__init__.py")

if(${BUILD_TESTING})
    add_subdirectory(tests)
endif()
