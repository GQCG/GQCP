list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/RHFElectricalResponseSolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vAP1roGElectricalResponseSolver.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
