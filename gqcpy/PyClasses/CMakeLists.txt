list(APPEND python_bindings_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Nucleus.cpp
        )

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
