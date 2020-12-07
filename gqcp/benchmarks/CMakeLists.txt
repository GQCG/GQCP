set(benchmark_target_sources)

add_subdirectory(ONVBasis)
add_subdirectory(QCMethod)


file(COPY data DESTINATION ${CMAKE_BINARY_DIR}/gqcp/benchmarks)  # make sure that the paths in the source files point to the correct data files


foreach(BENCHMARK_SOURCE ${benchmark_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(BENCHMARK_NAME ${BENCHMARK_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${BENCHMARK_NAME} ${BENCHMARK_SOURCE})

    # Configure (include headers and link libraries) the benchmark
    target_link_libraries(${BENCHMARK_NAME} PUBLIC gqcp)
    target_link_libraries(${BENCHMARK_NAME} PUBLIC benchmark::benchmark)  # link google benchmarks
endforeach()

# For program arguments, see https://github.com/google/benchmark
# --benchmark_counters_tabular=true (for nice format)
# --benchmark_out=<filename> (to save to file)
