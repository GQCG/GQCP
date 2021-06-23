list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GTOShell_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ScalarBasis_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ShellSet_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
