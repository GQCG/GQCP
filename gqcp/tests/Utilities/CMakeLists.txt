list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/miscellaneous_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/units_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
