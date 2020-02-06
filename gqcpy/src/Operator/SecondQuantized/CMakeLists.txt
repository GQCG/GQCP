list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/SQHamiltonian.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SQOneElectronOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SQTwoElectronOperator.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
