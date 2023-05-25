list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DiscreteDomain_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HubbardDomain_bindings.cpp
)

add_subdirectory(MullikenDomain)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
