cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(gqcp VERSION 0.2.0 LANGUAGES CXX)
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)  # make CMake look into the ./cmake/ folder for configuration files

# Supported options
# -----------------
option(BUILD_TESTS "Build tests" OFF)
option(BUILD_DOCS "Build the documentation using Doxygen" OFF)
option(BUILD_BENCHMARKS "Build benchmarks executables" OFF)
option(BUILD_DRIVERS "Build standard drivers" OFF)
option(EIGEN_USE_MKL "Eigen uses MKL" OFF)
option(BUILD_PYTHON_BINDINGS "Build the Python bindings" OFF)


# Setup libraries
# ---------------
add_library(gqcp SHARED "")

#if (BUILD_PYTHON_BINDINGS)
#    add_library(pygqcp SHARED "")
#endif()

# Look for supporting libraries
# -----------------------------
find_package(Git REQUIRED)
find_package(Eigen3 3.3.4 REQUIRED)
find_package(Int2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)
find_package(Spectra REQUIRED)
find_package(Cint REQUIRED NO_MODULE)
#if (BUILD_PYTHON_BINDINGS)
#    find_package(PythonInterp REQUIRED)
#    find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
#endif()
if(BUILD_PYTHON_BINDINGS)
    find_package(pybind11)
endif()

if(EIGEN_USE_MKL)
    set(BLA_VENDOR Intel10_64lp)
    find_package(BLAS REQUIRED)
endif()

if(BUILD_TESTS)
    #if (BUILD_PYTHON_BINDINGS)
    #    find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} numpy${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
    #else()
        find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)
    #endif()
else()
    #if (BUILD_PYTHON_BINDINGS)
    #    find_package(Boost REQUIRED COMPONENTS program_options python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} numpy${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})
    #else()
        find_package(Boost REQUIRED COMPONENTS program_options)
    #endif()
endif()

if (BUILD_DOCS)
    find_package(Doxygen REQUIRED dot)
endif()

if (BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
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
            Cint::Cint
            Spectra::Spectra
        )
target_include_directories(gqcp PUBLIC $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>)
target_link_libraries(gqcp PUBLIC ${Boost_LIBRARIES})

if (EIGEN_USE_MKL)
    target_include_directories(gqcp PUBLIC $<BUILD_INTERFACE:${BLAS_INCLUDE_DIR}>)
    target_link_libraries(gqcp PUBLIC ${BLAS_LIBRARIES})
    target_compile_options(gqcp PUBLIC -DEIGEN_USE_MKL_ALL -DMKL_LP64)
endif()

#if (BUILD_PYTHON_BINDINGS)
#    target_include_directories(pygqcp PUBLIC $<BUILD_INTERFACE:${PYTHON_INCLUDE_DIRS}>)
#    target_include_directories(pygqcp PUBLIC $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>)
#    target_link_libraries(pygqcp PUBLIC ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gqcp)
#    set_target_properties(pygqcp PROPERTIES SUFFIX ".so")
#    set_target_properties(pygqcp PROPERTIES PREFIX "")
#endif()


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

install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
        DESTINATION cmake
)


# Git related information
# -----------------------

# Get the long git SHA1 (https://stackoverflow.com/a/21028226/7930415)
execute_process(COMMAND "${GIT_EXECUTABLE}" describe --always --abbrev=40 --dirty
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_SHA1
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
