list(APPEND benchmark_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DOCI_CO_benchmark.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FCI_H-chain_benchmark.cpp
)

set(benchmark_target_sources ${benchmark_target_sources} PARENT_SCOPE)
