list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GTransformation_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RTransformation_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UTransformation_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UTransformationComponent_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
