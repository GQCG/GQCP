list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DiscreteDomainPartition_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ONVPartition_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
