list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/CI.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CIEnvironment.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
