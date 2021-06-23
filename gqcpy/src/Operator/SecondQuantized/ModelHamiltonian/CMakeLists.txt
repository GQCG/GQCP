list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/AdjacencyMatrix_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HoppingMatrix_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HubbardHamiltonian_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
