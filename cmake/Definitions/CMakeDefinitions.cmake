# Define the CMake environment


# Uppercase and lowercase names
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)


# Specify the library build type
set(LIBRARY_NAME ${PROJECT_NAME})
set(LIBRARY_TYPE SHARED)  # dynamic/shared library
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Building ${LIBRARY_NAME} in ${CMAKE_BUILD_TYPE} mode")


include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Options.cmake)


# Include all headers and source files
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Headers.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Sources.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Tests.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Drivers.cmake)

if (BUILD_BENCHMARKS)
    include(${CMAKE_SOURCE_DIR}/cmake/Definitions/Benchmarks.cmake)
endif()


# Find all packages
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Find)  # let CMake know that we have supplied our own FindXXX.cmake files
include(${CMAKE_SOURCE_DIR}/cmake/Find/FindPackages.cmake)


# Provide a function that takes includes/links this library to it
function(configure_executable EXECUTABLE_NAME)

    target_include_directories(${EXECUTABLE_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})  # include this project
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${LIBRARY_NAME})

endfunction(configure_executable)


include(${CMAKE_SOURCE_DIR}/cmake/Definitions/InstallDirectories.cmake)
