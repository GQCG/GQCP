list(APPEND python_bindings_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/DOCINewtonOrbitalOptimizer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DOCIRHF.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FCI.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FukuiDysonAnalysis.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Hubbard.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/MullikenConstrainedFCI.cpp
        )

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
