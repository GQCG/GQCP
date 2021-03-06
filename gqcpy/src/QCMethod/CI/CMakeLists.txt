list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/CI_factory_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CI_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CIEnvironment_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
