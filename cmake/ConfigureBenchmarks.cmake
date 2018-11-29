# In this file, we will provide a function that takes the name of a benchmark executable, and includes/links the appropriate libraries to it


function(configure_benchmark_executable EXECUTABLE_NAME)

    configure_executable(${EXECUTABLE_NAME})

    # Link google benchmarks
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC "-lbenchmark")


endfunction(configure_benchmark_executable)
