list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/OneDM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved1DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved2DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TwoDM_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
