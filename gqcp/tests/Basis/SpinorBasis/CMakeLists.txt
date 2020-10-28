list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GSpinorBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalSpace_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RSpinOrbitalBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleSpinorBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USpinOrbitalBasis_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
