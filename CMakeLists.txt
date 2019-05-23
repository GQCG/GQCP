cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(gqcp VERSION 0.2.0 LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_library(gqcp SHARED "")

# Look for supporting libraries
# -----------------------------

find_package(Git REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)
find_package(Eigen3 3.3.4 REQUIRED)
find_package(Int2 REQUIRED)
find_package(Spectra REQUIRED)
find_package(Cint REQUIRED)

if(EIGEN_USE_MKL_ALL)
    set(BLA_VENDOR Intel10_64lp)
    find_package(BLAS REQUIRED)
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
target_compile_features(gqcp PUBLIC cxx_std_11)

option(BUILD_DOCS "Build the documentation using Doxygen" OFF)
option(BUILD_BENCHMARKS "Build benchmarks executables" OFF)
option(BUILD_DRIVERS "Build standard drivers" OFF)
option(EIGEN_USE_MKL_ALL "Eigen uses MKL" ON)

add_subdirectory(src)
add_subdirectory(include)
enable_testing()
add_subdirectory(tests)

target_include_directories(gqcp
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        )

target_link_libraries(gqcp
        PUBLIC
            Eigen3::Eigen
            Int2::Int2
            Cint::Cint
            Spectra::Spectra
        )
target_include_directories(gqcp PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(gqcp PUBLIC ${Boost_LIBRARIES})
if (EIGEN_USE_MKL_ALL)
    target_include_directories(gqcp PUBLIC ${BLAS_INCLUDE_DIR})
    target_link_libraries(gqcp PUBLIC ${BLAS_LIBRARIES})
endif()

# Parse additional options
# ------------------------

# Generate documentation
if (BUILD_DOCS)
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/cmake/Parse/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(docs ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            VERBATIM)
endif()

if (BUILD_BENCHMARKS)
    add_subdirectory(${PROJECT_BENCHMARKS_FOLDER})
endif()

if (BUILD_DRIVERS)
    add_subdirectory(${PROJECT_DRIVERS_FOLDER})
endif()

# Setup installation directories
# ------------------------------

if (NOT INSTALLATION_PREFIX)
    set(INSTALLATION_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
set(PROJECT_INSTALL_DIR ${INSTALLATION_PREFIX}/${PROJECT_NAME_LOWERCASE})
set(INCLUDE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/include)
set(CMAKE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/cmake)
set(LIBRARY_INSTALL_DIR ${PROJECT_INSTALL_DIR}/lib)
set(BIN_INSTALL_DIR ${PROJECT_INSTALL_DIR}/bin)

# Install the library
# -------------------

## To specify that the library target should also be exported, we add the EXPORT option. This is used in conjuction with the install(EXPORT) command below
#install(TARGETS gqcp
#        EXPORT gqcp
#        LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})
#
#
## Install the header files
#install(DIRECTORY ${PROJECT_INCLUDE_FOLDER}/ DESTINATION ${INCLUDE_INSTALL_DIR})
#
#
## Export the target library into a ${PROJECT_NAME}Targets.cmake file.
## The file gqcpConfig.cmake includes this file, to be able to use this library with a find_package(template X.Y.Z) call
#install(EXPORT ${LIBRARY_NAME}
#        DESTINATION ${CMAKE_INSTALL_DIR}
#        FILE ${PROJECT_NAME}Targets.cmake)
#
## Install Config.cmake and ConfigVersion.cmake
#install(FILES
#        ${CMAKE_SOURCE_DIR}/cmake/Parse/${PROJECT_NAME}Config.cmake
#        ${CMAKE_SOURCE_DIR}/cmake/Parse/${PROJECT_NAME}ConfigVersion.cmake
#        DESTINATION ${CMAKE_INSTALL_DIR})



# From RPATH to full installation path
# https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
#set(CMAKE_INSTALL_RPATH "${LIBRARY_INSTALL_DIR}")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)