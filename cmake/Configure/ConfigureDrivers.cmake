# Configure and install the drivers


foreach(DRIVER_SOURCE ${PROJECT_DRIVERS_SOURCE_FILES})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(DRIVER_NAME ${DRIVER_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${DRIVER_NAME} ${DRIVER_SOURCE})

    # Configure (include headers and link libraries) the driver
    configure_executable(${DRIVER_NAME})

    # Install the driver in a separate location
    install(TARGETS ${DRIVER_NAME} RUNTIME DESTINATION ${BIN_INSTALL_DIR})

endforeach()
