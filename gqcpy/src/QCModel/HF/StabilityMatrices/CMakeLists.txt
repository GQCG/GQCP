list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GHFStabilityMatrices_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHFStabilityMatrices_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UHFStabilityMatrices_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
