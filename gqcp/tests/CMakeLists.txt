set(test_target_sources)

add_subdirectory(Basis)
add_subdirectory(ONVBasis)
add_subdirectory(Mathematical)
add_subdirectory(Molecule)
add_subdirectory(Operator)
add_subdirectory(Processing)
add_subdirectory(QCMethod)
add_subdirectory(Utilities)


file(COPY data DESTINATION ${CMAKE_BINARY_DIR}/gqcp/tests)  # make sure that the paths in the source files point to the correct data files


foreach(TEST_SOURCE ${test_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    # Configure (include headers and link libraries) the test
    target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${TEST_NAME} PUBLIC gqcp Boost::unit_test_framework)
    target_compile_definitions(${TEST_NAME} PRIVATE BOOST_TEST_DYN_LINK)

    add_test(
        NAME ${TEST_NAME} COMMAND ${TEST_NAME}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gqcp/tests
    )
endforeach()
