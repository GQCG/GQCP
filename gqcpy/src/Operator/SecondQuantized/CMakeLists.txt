list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/SQHamiltonian_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SQOneElectronOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SQTwoElectronOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USQHamiltonian_bindings.cpp
)

add_subdirectory(ModelHamiltonian)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
