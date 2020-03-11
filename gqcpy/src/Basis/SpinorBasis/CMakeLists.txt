list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GSpinorBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RSpinorBasis.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
