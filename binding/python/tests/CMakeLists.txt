add_test(NAME python-tests
    COMMAND ${PYTHON_EXECUTABLE} -m pyTests # Or just COMMAND pytest
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)