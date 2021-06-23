list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/ERJacobiLocalizer_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ERNewtonLocalizer_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
