# Find all packages


find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Eigen3 3.3.4 REQUIRED)
find_package(Libint2 REQUIRED)
find_package(Spectra REQUIRED)

if (BUILD_DOCS)
    find_package(Doxygen REQUIRED dot)
endif()

if (USE_MKL)
    find_package(MKL REQUIRED)
endif()

if (BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
endif()
