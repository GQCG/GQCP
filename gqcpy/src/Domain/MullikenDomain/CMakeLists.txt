list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GMullikenDomain_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RMullikenDomain_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UMullikenDomainComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UMullikenDomain_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
