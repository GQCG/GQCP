list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/Simple1DM_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Simple2DM_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved1DM_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved2DM_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
