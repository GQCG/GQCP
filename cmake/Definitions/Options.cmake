# This CMake-file takes care of defining all the options


option(BUILD_DOCS "Build the documentation using Doxygen" OFF)
option(BUILD_BENCHMARKS "Build benchmarks executables" OFF)
option(BUILD_DRIVERS "Build standard drivers" OFF)

option(USE_MKL "Find and use MKL for BLAS" OFF)
