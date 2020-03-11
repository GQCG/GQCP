list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roG.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGGeminalCoefficients.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/vAP1roG.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
