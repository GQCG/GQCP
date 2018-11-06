# benchmark-level CMakeLists.txt

# For every *.cpp-file (every file corresponding to an executable) ...
foreach(BENCH_SOURCE ${PROJECT_BENCH_SOURCE_FILES})


    link_libraries("-lbenchmark")
    #add_compile_options("-benchmark_counters_tabular=true")

    # ... extract the filename without extension (NAME_WE) as a name for our executable ...
    get_filename_component(EXE_NAME ${BENCH_SOURCE} NAME_WE)

    # ... add an executable based on the source ...
    add_executable(${EXE_NAME} ${BENCH_SOURCE})

    # ... configure (include headers and link libraries) the executable ...
    configure_executable(${EXE_NAME})
endforeach()
