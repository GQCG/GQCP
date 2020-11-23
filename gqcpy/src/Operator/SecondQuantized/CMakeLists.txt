list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/SQHamiltonian_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RSQOneElectronOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RSQTwoElectronOperator_bindings.cpp
)

add_subdirectory(ModelHamiltonian)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
