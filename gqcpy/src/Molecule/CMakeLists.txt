list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/Molecule_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearFramework_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Nucleus_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
