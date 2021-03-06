list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/AdjacencyMatrix_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HoppingMatrix_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
