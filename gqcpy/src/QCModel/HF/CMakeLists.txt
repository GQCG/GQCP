list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GHF_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHF_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UHF_bindings.cpp
)

add_subdirectory(StabilityMatrices)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
