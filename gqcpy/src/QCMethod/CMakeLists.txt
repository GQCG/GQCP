add_subdirectory(CC)
add_subdirectory(CI)
add_subdirectory(Geminals)
add_subdirectory(HF)
add_subdirectory(NOCI)

list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/QCStructure_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
