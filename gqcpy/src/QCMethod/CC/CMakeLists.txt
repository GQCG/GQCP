list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/CCSD_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CCSDEnvironment_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CCSDSolver_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
