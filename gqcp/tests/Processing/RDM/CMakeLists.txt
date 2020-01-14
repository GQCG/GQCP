list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DOCIRDMBuilder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FCIRDMBuilder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FrozenCoreDOCIRDMBuilder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FrozenCoreFCIRDMBuilder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OneRDM_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RDMCalculator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SelectedRDMBuilder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinUnresolvedFCIRDMBuilder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TwoRDM_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
