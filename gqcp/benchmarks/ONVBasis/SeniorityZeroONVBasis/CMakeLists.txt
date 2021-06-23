list(APPEND benchmark_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/SeniorityZeroONVBasis_RSQHamiltonian_dense_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SeniorityZeroONVBasis_RSQHamiltonian_matvec_benchmark.cpp
)

set(benchmark_target_sources ${benchmark_target_sources} PARENT_SCOPE)
