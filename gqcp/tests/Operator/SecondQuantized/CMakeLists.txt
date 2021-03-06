add_subdirectory(ModelHamiltonian)

list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/EvaluatableScalarRSQOneElectronOperator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleSQOneElectronOperator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SQHamiltonian_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleSQTwoElectronOperator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USQOneElectronOperator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USQTwoElectronOperator_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
