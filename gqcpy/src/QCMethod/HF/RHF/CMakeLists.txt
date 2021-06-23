list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DiagonalRHFFockMatrixObjective_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHF_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHFSCFEnvironment_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHFSCFSolver_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
