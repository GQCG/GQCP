list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/atomic_decomposition_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/expectation_values_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/properties_test.cpp
        )
set(test_target_sources ${test_target_sources} PARENT_SCOPE)
