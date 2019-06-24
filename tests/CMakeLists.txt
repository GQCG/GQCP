set(test_target_sources)

list(APPEND test_target_sources
        Atom_test.cpp
        elements_test.cpp
        HoppingMatrix_test.cpp
        Molecule_test.cpp
        RMP2_test.cpp
        units_test.cpp
        )

add_subdirectory(Basis)
add_subdirectory(CISolver)
add_subdirectory(FockSpace)
add_subdirectory(Geminals)
add_subdirectory(HamiltonianBuilder)
add_subdirectory(HamiltonianParameters)
add_subdirectory(Mathematical)
add_subdirectory(Operator)
add_subdirectory(OrbitalOptimization)
add_subdirectory(Properties)
add_subdirectory(RDM)
add_subdirectory(RHF)
add_subdirectory(utilities)
add_subdirectory(WaveFunction)

file(COPY data DESTINATION ${CMAKE_BINARY_DIR})  # copy the necessary data files to the out-of-source build directory

foreach(TEST_SOURCE ${test_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${TEST_NAME} ${TEST_SOURCE})

    # Configure (include headers and link libraries) the test
    target_include_directories(${TEST_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${TEST_NAME} PUBLIC gqcp)
    target_compile_definitions(${TEST_NAME} PRIVATE BOOST_TEST_DYN_LINK)

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR})  # the working directory is the out-of-source build directory

endforeach()
