#
# Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

add_test(
  NAME python-tests
  COMMAND ${PYTHON_EXECUTABLE} -m pyTests # Or just COMMAND pytest
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
