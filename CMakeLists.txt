# top-level CMakeLists.txt

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(gqcp
        VERSION 0.1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)


# Establish the CMake environment - define all CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/Definitions/CMakeDefinitions.cmake)


# Parse all *.in-files
include(${CMAKE_SOURCE_DIR}/cmake/Parse/ParseInputFiles.cmake)


# Configure and build
add_subdirectory(${PROJECT_SOURCE_FOLDER})
add_subdirectory(${PROJECT_DRIVERS_FOLDER})

if (BUILD_BENCHMARKS)
    add_subdirectory(${PROJECT_BENCHMARKS_FOLDER})
endif()

enable_testing()
add_subdirectory(${PROJECT_TESTS_FOLDER})
