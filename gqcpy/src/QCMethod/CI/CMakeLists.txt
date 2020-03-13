list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/CI_factory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CI.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CIEnvironment.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DOCINewtonOrbitalOptimizer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DOCINewtonOrbitalOptimizer_factory.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
