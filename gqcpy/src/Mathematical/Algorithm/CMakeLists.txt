list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/Algorithm.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionalStep.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IterativeAlgorithm.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
