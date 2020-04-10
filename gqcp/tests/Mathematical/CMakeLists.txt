list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/ScalarFunction_test.cpp
)

add_subdirectory(Grid)
add_subdirectory(Optimization)
add_subdirectory(Representation)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
