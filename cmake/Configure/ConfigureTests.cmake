# Configure the test executables


file(COPY data DESTINATION ${CMAKE_BINARY_DIR})  # copy the necessary data files to the out-of-source build directory



foreach(TEST_SOURCE ${PROJECT_TEST_SOURCE_FILES})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    # Configure (include headers and link libraries) the test
    configure_executable(${TEST_NAME})
    target_include_directories(${TEST_NAME} PUBLIC ${Spectra_INCLUDE_DIRS})

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR})  # the working directory is the out-of-source build directory

endforeach()
