list(APPEND test_target_sources
    # ${CMAKE_CURRENT_SOURCE_DIR}/FrozenCoreDOCIRDMBuilder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GeneralDMCalculator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SeniorityZeroDMCalculator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedDMCalculator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedFrozenDMCalculator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolvedSelectedDMCalculator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinUnresolvedDMCalculator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TwoRDM_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
