list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/MixedUSQTwoElectronOperatorComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PureUSQTwoElectronOperatorComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RSQOneElectronOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RSQTwoElectronOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SQHamiltonian_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USQOneElectronOperator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USQOneElectronOperatorComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USQTwoElectronOperator_bindings.cpp
)

add_subdirectory(ModelHamiltonian)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
