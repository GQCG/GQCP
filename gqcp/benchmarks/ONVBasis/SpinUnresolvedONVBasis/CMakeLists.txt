list(APPEND benchmark_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinUnresolvedONVBasis_GSQOneElectronOperator_matvec_benchmark.cpp
)

set(benchmark_target_sources ${benchmark_target_sources} PARENT_SCOPE)
