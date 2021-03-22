list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearRepulsionOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearDipoleOperator_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
