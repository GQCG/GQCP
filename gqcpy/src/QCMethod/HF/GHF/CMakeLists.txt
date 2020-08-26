list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GHF_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GHFSCFEnvironment_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GHFSCFSolver_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)