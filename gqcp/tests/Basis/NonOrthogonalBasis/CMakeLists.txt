list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GNonOrthogonalStateBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RNonOrthogonalStateBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UNonOrthogonalStateBasis_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
