# The entry-point for the CMake configuration of the GQCP project
# The targets in the GQCP project are:
#   - the gqcp C++ library
#   - the executables for the unit tests of the C++ library (optional)
#   - the executables for the benchmarking of the C++ library (optional)
#   - the Doxygen-generated documentation of the C++ library (optional)
#   - the gqcpy Python bindings to the C++ library (optional)

cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(GQCP VERSION 0.1.0 LANGUAGES CXX)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)  # make CMake look into the ./cmake/ folder for configuration files

# Look for the dependencies
find_package(Git REQUIRED)
find_package(Eigen3 3.3.4 REQUIRED)
find_package(Int2 REQUIRED)
find_package(cint REQUIRED)
find_package(MKL REQUIRED)


# Get the latest commit hash
execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_SHA1
    OUTPUT_STRIP_TRAILING_WHITESPACE
)


# Add the versioning as a separate library. Whenever the version changes, only the version.cpp source file will be recompiled
# and only gqcp needs to be relinked. This way, gqcp does not need to be recompiled if the version should be modified.
# Add configuration file containing version of current library.
configure_file(${CMAKE_SOURCE_DIR}/cmake/version.hpp.in gqcp/version/version.hpp)
configure_file(${CMAKE_SOURCE_DIR}/cmake/version.cpp.in gqcp/version/version.cpp)

add_library(gqcp_version STATIC ${CMAKE_CURRENT_BINARY_DIR}/gqcp/version/version.cpp)
target_compile_features(gqcp_version PUBLIC cxx_std_14)
target_include_directories(gqcp_version
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/gqcp/version>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/gqcp/version>
        $<INSTALL_INTERFACE:${CMAKE_PREFIX_PATH}/include/gqcp/version>
)

# List the supported options.
option(BUILD_TESTS "Build the C++ library tests" OFF)
option(BUILD_BENCHMARKS "Build the executables used for benchmarking the C++ library" OFF)
option(BUILD_DOCS "Build the documentation of the C++ library using Doxygen" OFF)
option(BUILD_PYTHON_BINDINGS "Build the Python bindings for the C++ library" OFF)
option(OPTIMIZE_FOR_NATIVE "Build with -march=native" OFF)


if (BUILD_TESTS)
    find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)
else()
    find_package(Boost REQUIRED COMPONENTS program_options)
endif()

if (BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
endif()

if (BUILD_DOCS)
    find_package(Doxygen REQUIRED dot)
endif()

# Look for the dependencies: find Python
# Pybind11 conflicts with the CMake Python detection (https://pybind11.readthedocs.io/en/stable/faq.html#inconsistent-detection-of-python-version-in-cmake-and-pybind11)
# If the Python bindings are built, use Pybind11 to detect Python, else let standard CMake do it
if (BUILD_PYTHON_BINDINGS)
    # Passing PYTHON_EXECUTABLE as a variable into CMake helps PythonInterp (which is called through PyBind11) find the correct Python executable and libraries
    find_package(pybind11 REQUIRED)
    set(GQCPYTHON_INTERPRETER ${PYTHON_EXECUTABLE})  # PyBind11 uses find_package(PythonInterp), which sets this variable
else()
    find_package(Python3 COMPONENTS Interpreter)
    set(GQCPYTHON_INTERPRETER ${Python3_EXECUTABLE})
endif()


# Setup the C++ library

# Improve the handling of RPATH settings (https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(gqcp)



# Parse the options: tests for the C++ library
if (BUILD_TESTS)
    enable_testing()
    add_subdirectory(gqcp/tests)
endif()


# Parse the options: executables for the benchmarks of the C++ library
if (BUILD_BENCHMARKS)
    add_subdirectory(gqcp/benchmarks)
endif()


# Parse the options: documentation of the C++ library
if (BUILD_DOCS)
    set(DOXYGEN_IN cmake/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(
        docs ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        VERBATIM
    )
endif()


# Parse the options: Python bindings
if (BUILD_PYTHON_BINDINGS)
    add_subdirectory(gqcpy)
endif()



# Install the C++ library (relative to the CMAKE_INSTALL_PREFIX)
include(GNUInstallDirs)

install(
    TARGETS gqcp_version
    EXPORT gqcp-targets
    LIBRARY DESTINATION lib/gqcp
    ARCHIVE DESTINATION lib/gqcp
    RUNTIME DESTINATION bin/gqcp
    INCLUDES DESTINATION include/gqcp/version
)

install(DIRECTORY ${CMAKE_BINARY_DIR}/gqcp/version/ DESTINATION include/gqcp/version)

install(
    TARGETS gqcp
    EXPORT gqcp-targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include/gqcp
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/gqcp/include/ DESTINATION include/gqcp)

install(
    EXPORT gqcp-targets
    FILE gqcp-targets.cmake
    NAMESPACE GQCP::
    DESTINATION cmake
)


# Configure the GQCP package for downstream CMake importing
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/gqcp-config.cmake.in
    ${CMAKE_BINARY_DIR}/cmake/gqcp-config.cmake
    INSTALL_DESTINATION cmake
)

write_basic_package_version_file(
    ${CMAKE_BINARY_DIR}/cmake/gqcp-config-version.cmake
    VERSION ${gqcp_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES
        ${CMAKE_BINARY_DIR}/cmake/gqcp-config.cmake
        ${CMAKE_BINARY_DIR}/cmake/gqcp-config-version.cmake
        ${CMAKE_SOURCE_DIR}/cmake/FindEigen3.cmake
        ${CMAKE_SOURCE_DIR}/cmake/FindInt2.cmake
        ${CMAKE_SOURCE_DIR}/cmake/FindMKL.cmake
    DESTINATION cmake
)


# Install the Python bindings
if (BUILD_PYTHON_BINDINGS)
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} \"\${CMAKE_BINARY_DIR}/gqcpy/setup.py\" install --prefix ${CMAKE_INSTALL_PREFIX} --force)")
endif()
