list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roG_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGLagrangianNewtonOrbitalOptimizer_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PSEnvironment_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vAP1roG_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
