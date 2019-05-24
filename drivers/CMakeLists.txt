set(driver_target_sources)

list(APPEND driver_target_sources
        fci_lowdin.cpp
        hubbard.cpp
        oo_doci.cpp
        )

foreach(DRIVER_SOURCE ${driver_target_sources})
    # Extract the filename without extension (NAME_WE) as a name for our executable
    get_filename_component(DRIVER_NAME ${DRIVER_SOURCE} NAME_WE)

    # Add an executable based on the source
    add_executable(${DRIVER_NAME} ${DRIVER_SOURCE})

    # Configure (include headers and link libraries) the driver
    target_include_directories(${DRIVER_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${DRIVER_NAME} PUBLIC gqcp)

    # Install the driver in a separate location
    install(TARGETS ${DRIVER_NAME} RUNTIME DESTINATION bin)

endforeach()