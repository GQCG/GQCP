list(APPEND benchmark_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedONVBasis_HubbardHamiltonian_dense_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedONVBasis_HubbardHamiltonian_matvec_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedONVBasis_RSQHamiltonian_dense_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedONVBasis_RSQHamiltonian_matvec_benchmark.cpp
)

set(benchmark_target_sources ${benchmark_target_sources} PARENT_SCOPE)
