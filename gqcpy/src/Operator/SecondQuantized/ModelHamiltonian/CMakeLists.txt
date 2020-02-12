list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/HoppingMatrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HubbardHamiltonian.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
