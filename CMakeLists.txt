# src-level CMakeLists.txt

cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

# Specify the project, its version and the used language
project(gqcp
        VERSION 0.1.0
        LANGUAGES CXX)

# Option to build the documentation
option(BUILD_DOCS "Build the documentation using Doxygen" OFF)

# Option to build benchmarks
option(BUILD_BENCHMARKS "Build the benchmarks" OFF)

# To enable the C++11 standard, we have to add flags like -std=gnu+11 or so. CMake can take care of this:
set(CMAKE_CXX_STANDARD 11)

# Set some compiler optimizations (https://stackoverflow.com/a/41361741/7930415)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if (“${CMAKE_CXX_COMPILER_ID}” STREQUAL “Clang”)
    set(CMAKE_CXX_FLAGS_DEBUG "-g -march=native -m64 -pipe")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -m64 -pipe")
elseif (“${CMAKE_CXX_COMPILER_ID}” STREQUAL “GNU”)
    set(CMAKE_CXX_FLAGS_DEBUG "-g -march=native -m64 -pipe")  
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -m64 -pipe")
  elseif (“${CMAKE_CXX_COMPILER_ID}” STREQUAL “Intel”)
    set(CMAKE_CXX_FLAGS_RELEASE "-g -xHost -m64 -pipe")    
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -xHost -m64 -pipe")
endif()

# Set CMake variables
include(${CMAKE_SOURCE_DIR}/cmake/SetCMakeEnvironment.cmake)

message(STATUS "Building ${LIBRARY_NAME} in ${CMAKE_BUILD_TYPE} mode")

# Find all required packages. Update CMake's MODULE_PATH to let CMake know that we have supplied our own FindXXX.cmake files in the cmake directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include(${CMAKE_SOURCE_DIR}/cmake/FindPackages.cmake)

# Add the sources for the library
add_subdirectory(${PROJECT_SOURCE_FOLDER})

# Add the sources for the executables
add_subdirectory(${PROJECT_EXECUTABLES_FOLDER})

# Add the sources for the benchmarks
if(BUILD_BENCHMARKS)
    add_subdirectory(${PROJECT_BENCHMARKS_FOLDER})
endif(BUILD_BENCHMARKS)

# Set up testing
enable_testing()
add_subdirectory(${PROJECT_TESTS_FOLDER})
