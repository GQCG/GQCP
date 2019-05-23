# Find all packages
find_package(Git REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options unit_test_framework)
find_package(Eigen3 3.3.4 REQUIRED)
find_package(Int2 REQUIRED)
find_package(Cint REQUIRED)

if (BUILD_DOCS)
    find_package(Doxygen REQUIRED dot)
endif()

if(EIGEN_USE_MKL_ALL)
    set(BLA_VENDOR Intel10_64lp)
    find_package(BLAS REQUIRED)
endif()

if (BUILD_BENCHMARKS)
    find_package(benchmark REQUIRED)
endif()
