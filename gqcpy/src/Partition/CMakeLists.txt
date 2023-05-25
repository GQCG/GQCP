list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DiscreteDomainPartition_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ONVPartition_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedElectronPartition_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinUnresolvedElectronPartition_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
