list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GLowdinPairingBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RLowdinPairingBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ULowdinPairingBasis_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
