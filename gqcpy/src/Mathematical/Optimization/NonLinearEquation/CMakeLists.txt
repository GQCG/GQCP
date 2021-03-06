list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/NonLinearEquationEnvironment_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NonLinearEquationSolver_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
