# In this file, we will provide a function that takes the name of a benchmark executable, and includes/links the appropriate libraries to it


function(configure_benchmark_executable EXECUTABLE_NAME)

    # Include this project
    target_include_directories(${EXECUTABLE_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${LIBRARY_NAME})

    # Link google benchmarks
    link_libraries("-lbenchmark")

endfunction(configure_benchmark_executable)
