list(APPEND python_bindings_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/bindings.cpp
    )

add_subdirectory(QCMethod)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
