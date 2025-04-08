set(sandbox_target_sources)

list(APPEND sandbox_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/test_sandbox.cpp
)

set(sandbox_target_sources ${sandbox_target_sources} PARENT_SCOPE)

# file(COPY data DESTINATION ${CMAKE_BINARY_DIR}/gqcp/sandbox)  # make sure that the paths in the source files point to the correct data files

foreach(SANDBOX_SOURCE ${sandbox_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(SANDBOX_NAME ${SANDBOX_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${SANDBOX_NAME} ${SANDBOX_SOURCE})

    # Configure (include headers and link libraries) the benchmark
    target_link_libraries(${SANDBOX_NAME} PUBLIC gqcp)
endforeach()
