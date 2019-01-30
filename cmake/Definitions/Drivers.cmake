# List all drivers


# Find the drivers folder
set(PROJECT_DRIVERS_FOLDER ${CMAKE_SOURCE_DIR}/drivers)


# Find the source files for the drivers
set(PROJECT_DRIVERS_SOURCE_FILES
        ${PROJECT_DRIVERS_FOLDER}/fci_lowdin.cpp
        ${PROJECT_DRIVERS_FOLDER}/hubbard.cpp
        ${PROJECT_DRIVERS_FOLDER}/oo_doci.cpp
    )
