list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/CoulombRepulsionOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ElectronicDipoleOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/KineticOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearAttractionOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearRepulsionOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearDipoleOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OverlapOperator_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
