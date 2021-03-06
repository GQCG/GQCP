list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DysonOrbital_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHFElectricalResponseSolver_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vAP1roGElectricalResponseSolver_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
