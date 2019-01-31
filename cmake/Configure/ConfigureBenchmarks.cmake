# Provide a function takes the name of a target benchmark executable, and includes/links the appropriate libraries to it


foreach(BENCHMARK_SOURCE ${PROJECT_BENCHMARK_SOURCE_FILES})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(BENCHMARK_NAME ${BENCHMARK_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${BENCHMARK_NAME} ${BENCHMARK_SOURCE})

    # Configure (include headers and link libraries) the benchmark ...
    configure_executable(${BENCHMARK_NAME})
    target_link_libraries(${BENCHMARK_NAME} PUBLIC "-lbenchmark")  # link google benchmarks

    # Install the benchmark in a separate location
    install(TARGETS ${BENCHMARK_NAME} RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endforeach()


# For program arguments, see https://github.com/google/benchmark
# --benchmark_counters_tabular=true (for nice format)
# --benchmark_out=<filename> (to save to file)
