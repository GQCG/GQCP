list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/CartesianDirection_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CartesianExponents_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CartesianGTO_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
