list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/interfaces.hpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
