list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/SQHamiltonian_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SQOneElectronOperator_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SQTwoElectronOperator_test.cpp
        )

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
