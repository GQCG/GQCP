list(APPEND python_bindings_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/DOCINewtonOrbitalOptimizer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FCI.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Hubbard.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MullikenConstrainedFCI.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Nucleus.cpp
        )

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
