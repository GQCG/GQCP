list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DiagonalRHFFockMatrixObjective.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHFSCFEnvironment.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHFSCFSolver.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
