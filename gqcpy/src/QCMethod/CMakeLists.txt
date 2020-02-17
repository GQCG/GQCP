add_subdirectory(CI)
add_subdirectory(HF)

list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/QCStructure.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
