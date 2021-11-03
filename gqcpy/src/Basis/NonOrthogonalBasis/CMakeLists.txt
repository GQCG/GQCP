list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GNonOrthogonalStateBasis_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RNonOrthogonalStateBasis_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UNonOrthogonalStateBasis_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
