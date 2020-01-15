# Entry point for the CMake configuration of the gqcpy Python bindings

# Add the sources for the Python bindings
set(python_bindings_sources)

list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/gqcpy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp
)

add_subdirectory(src)


# Configure the Python bindings
pybind11_add_module(gqcpy MODULE ${python_bindings_sources})

target_link_libraries(gqcpy PUBLIC gqcp)
set_target_properties(gqcpy
    PROPERTIES
        SUFFIX ".so"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gqcpy"
)
configure_file(__init__.py.in ${CMAKE_BINARY_DIR}/gqcpy/__init__.py)
configure_file(setup.py.in ${CMAKE_BINARY_DIR}/gqcpy/setup.py)


# Configure the Python tests
set(python_test_target_sources)
add_subdirectory(tests)

file(COPY tests/data DESTINATION ${CMAKE_BINARY_DIR}/gqcpy/tests)  # copy the necessary data files to the out-of-source build directory

foreach(TEST_SOURCE ${python_test_target_sources})
    # Extract the filename
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME)

    # Copy the file to the out-of-source build gqcpy tests directory
    configure_file(${TEST_SOURCE} ${CMAKE_BINARY_DIR}/gqcpy/tests/${TEST_NAME} COPYONLY)

    # Add test as python run command
    add_test(NAME ${TEST_NAME} COMMAND ${GQCPYTHON_INTERPRETER} -m unittest ${CMAKE_BINARY_DIR}/gqcpy/tests/${TEST_NAME}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gqcpy/tests)
endforeach()
