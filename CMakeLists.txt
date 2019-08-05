cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(gqcp VERSION 0.2.0 LANGUAGES CXX)
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)

# Improved handling of RPATH settings
# https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)  # make CMake look into the ./cmake/ folder for configuration files

add_library(gqcp SHARED "")

# Supported options
# -----------------
option(BUILD_TESTS "Build tests" OFF)
option(BUILD_DOCS "Build the documentation using Doxygen" OFF)
option(BUILD_BENCHMARKS "Build benchmarks executables" OFF)
option(BUILD_DRIVERS "Build standard drivers" OFF)
option(BUILD_PYTHON_BINDINGS "Build the Python bindings" OFF)

# Look for supporting libraries
# -----------------------------
find_package(Git REQUIRED)
find_package(Eigen3 3.3.4 REQUIRED)
find_package(Int2 REQUIRED)
find_package(Spectra REQUIRED)
find_package(cint REQUIRED)
find_package(MKL REQUIRED)
find_package(Python3 COMPONENTS Interpreter)

if(BUILD_TESTS)
    find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)
else()
    find_package(Boost REQUIRED COMPONENTS program_options)
endif()

if (BUILD_DOCS)
    find_package(Doxygen REQUIRED dot)
endif()

if (BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
endif()

if(BUILD_PYTHON_BINDINGS)
    find_package(PythonInterp 3 REQUIRED)
    find_package(pybind11 REQUIRED)
    set(python_bindings_sources)
endif()


# Setup library
# -------------
target_compile_options(gqcp PUBLIC -m64 -pipe)
target_compile_options(gqcp PUBLIC "$<IF:$<CONFIG:Debug>,-g,-O2>")
target_compile_options(gqcp PUBLIC "$<IF:$<CXX_COMPILER_ID:Intel>,-xHost,-march=native>")
target_compile_options(gqcp PUBLIC "$<$<CXX_COMPILER_ID:GNU>:-pthread>")
target_link_options(gqcp PUBLIC "$<$<CXX_COMPILER_ID:GNU>:-pthread>")
target_compile_features(gqcp PUBLIC cxx_std_11)

add_subdirectory(src)
add_subdirectory(include)

# Get the long git SHA1 (https://stackoverflow.com/a/21028226/7930415)
execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always --abbrev=40 --dirty
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_SHA1
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file(cmake/version.hpp.in include/version.hpp)

target_include_directories(gqcp
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
        PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        )

target_link_libraries(gqcp
        PUBLIC
            Eigen3::Eigen
            Int2::Int2
            cint::cint
            Spectra::Spectra
            MKL::MKL
        )
target_include_directories(gqcp PUBLIC $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>)
target_link_libraries(gqcp PUBLIC ${Boost_LIBRARIES})
target_compile_options(gqcp PUBLIC -DEIGEN_USE_MKL_ALL -DMKL_LP64)


# Setup Python Bindings
# ---------------------

if (BUILD_PYTHON_BINDINGS)
    add_subdirectory(gqcpy)

    pybind11_add_module(gqcpy MODULE ${python_bindings_sources})
    
    target_link_libraries(gqcpy PUBLIC gqcp)
    set_target_properties(gqcpy PROPERTIES
            SUFFIX ".so"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gqcpy")
    configure_file(gqcpy/__init__.py.in ${CMAKE_BINARY_DIR}/gqcpy/__init__.py)
    configure_file(gqcpy/setup.py.in ${CMAKE_BINARY_DIR}/gqcpy/setup.py)
endif()


# Parse additional options
# ------------------------
if (BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()


# Generate documentation
# ----------------------
if (BUILD_DOCS)
    set(DOXYGEN_IN cmake/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(docs ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            VERBATIM)
endif()

if (BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if (BUILD_DRIVERS)
    add_subdirectory(drivers)
endif()


# Install the library (relative to the CMAKE_INSTALL_PREFIX)
# ----------------------------------------------------------

include(GNUInstallDirs)

install(TARGETS gqcp
        EXPORT gqcp-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/ DESTINATION include)
install(FILES ${CMAKE_BINARY_DIR}/include/version.hpp DESTINATION include)

install(EXPORT gqcp-targets
        FILE gqcp-targets.cmake
        NAMESPACE gqcp::
        DESTINATION cmake
        )

# Configure package
# -----------------

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
        DESTINATION cmake
)


# Install Python bindings
# -----------------------

if(BUILD_PYTHON_BINDINGS)
    install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} \"\${CMAKE_BINARY_DIR}/gqcpy/setup.py\" install --prefix ${CMAKE_INSTALL_PREFIX} --force)")
endif()
