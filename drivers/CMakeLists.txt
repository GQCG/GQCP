# drivers-level CMakeLists.txt

foreach(DRIVER_SOURCE ${PROJECT_DRIVER_SOURCE_FILES})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(DRIVER_NAME ${DRIVER_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${DRIVER_NAME} ${DRIVER_SOURCE})

    # Configure (include headers and link libraries) the executable
    configure_driver(${DRIVER_NAME})
endforeach()
