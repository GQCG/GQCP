list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/NonLinearEquationEnvironment.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NonLinearEquationSolver.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
