list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/NOCI_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NOCI_factory_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NOCIEnvironment_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
