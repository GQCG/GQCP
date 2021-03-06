list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/CubicGrid_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Field_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WeightedGrid_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
