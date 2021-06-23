list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearDipoleOperator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearRepulsionOperator_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
