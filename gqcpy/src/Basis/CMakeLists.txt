add_subdirectory(SpinorBasis)

list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/transform_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
