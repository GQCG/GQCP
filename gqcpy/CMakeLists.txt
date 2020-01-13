# Add the sources for the Python bindings
list(APPEND python_bindings_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp
    )

add_subdirectory(Basis)
add_subdirectory(Molecule)
add_subdirectory(Operator)
add_subdirectory(QCMethod)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)


# Configure the tests
set(python_test_target_sources)
add_subdirectory(tests)

file(COPY tests/data DESTINATION ${CMAKE_BINARY_DIR}/gqcpy/tests)  # copy the necessary data files to the out-of-source build directory

foreach(TEST_SOURCE ${python_test_target_sources})
    # Extract the filename 
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME)

    # Copy the file to the out-of-source build gqcpy tests directory
    configure_file(${TEST_SOURCE} ${CMAKE_BINARY_DIR}/gqcpy/tests/${TEST_NAME} COPYONLY)

    # Add test as python run command
    add_test(NAME ${TEST_NAME} COMMAND Python3::Interpreter -m unittest ${CMAKE_BINARY_DIR}/gqcpy/tests/${TEST_NAME}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gqcpy/tests)
endforeach()
