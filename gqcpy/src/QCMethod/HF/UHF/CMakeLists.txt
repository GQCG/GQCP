list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/UHF_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UHFSCFEnvironment_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UHFSCFSolver_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
