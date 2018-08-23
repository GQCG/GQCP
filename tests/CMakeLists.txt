# test-level CMakeLists.txt

# For every *.cpp-file (every file corresponding to one set of tests) ...
foreach(TEST_SOURCE ${PROJECT_TEST_SOURCE_FILES})
    # ... extract the filename without extension (NAME_WE) as a name for our executable ...
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    # ... add an executable based on the test source ...
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    # ... configure (include headers and link libraries) the test ...
    configure_executable(${TEST_NAME})

    # ... and finally make a test
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})  # the working directory will be the build directory
endforeach()
