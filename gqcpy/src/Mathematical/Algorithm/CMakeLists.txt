list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/Algorithm_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionalStep_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IterativeAlgorithm_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
