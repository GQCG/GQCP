list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GLowdinPairingBasis_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RLowdinPairingBasis_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ULowdinPairingBasis_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
