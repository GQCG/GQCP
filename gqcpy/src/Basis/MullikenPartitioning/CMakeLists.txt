list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/RMullikenPartitioning_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UMullikenPartitioningComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UMullikenPartitioning_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
