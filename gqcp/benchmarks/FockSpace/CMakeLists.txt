list(APPEND benchmark_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/fockspace_one_electron_matvec.cpp
)

set(benchmark_target_sources ${benchmark_target_sources} PARENT_SCOPE)
