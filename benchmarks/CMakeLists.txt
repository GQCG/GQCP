
set(benchmark_target_sources)

add_subdirectory(DOCI)
add_subdirectory(FCI)
add_subdirectory(Hubbard)

# Provide a function takes the name of a target benchmark executable, and includes/links the appropriate libraries to it

foreach(BENCHMARK_SOURCE ${benchmark_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(BENCHMARK_NAME ${BENCHMARK_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${BENCHMARK_NAME} ${BENCHMARK_SOURCE})

    # Configure (include headers and link libraries) the benchmark ...
    target_include_directories(${BENCHMARK_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${BENCHMARK_NAME} PUBLIC gqcp)
    target_link_libraries(${BENCHMARK_NAME} PUBLIC "-lbenchmark")  # link google benchmarks

    # Install the benchmark in a separate location
    #install(TARGETS ${BENCHMARK_NAME} RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endforeach()


# For program arguments, see https://github.com/google/benchmark
# --benchmark_counters_tabular=true (for nice format)
# --benchmark_out=<filename> (to save to file)
