list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/EigenproblemEnvironment.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EigenproblemSolver.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
