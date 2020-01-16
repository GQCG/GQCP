list(APPEND benchmark_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/fci_hchain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fci_matrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fci_matvec.cpp
)

set(benchmark_target_sources ${benchmark_target_sources} PARENT_SCOPE)
