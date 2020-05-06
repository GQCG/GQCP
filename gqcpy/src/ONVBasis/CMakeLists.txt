list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedONVBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedONV.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SeniorityZeroONVBasis.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
