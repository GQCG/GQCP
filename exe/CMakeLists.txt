# exe-level CMakeLists.txt


foreach(EXE_SOURCE ${PROJECT_EXE_SOURCE_FILES})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(EXE_NAME ${EXE_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${EXE_NAME} ${EXE_SOURCE})

    # Configure (include headers and link libraries) the executable
    configure_driver(${EXE_NAME})
endforeach()
