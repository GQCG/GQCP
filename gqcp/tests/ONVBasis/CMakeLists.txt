list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/ONVBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FrozenONVBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FrozenProductONVBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ONV_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ProductONVBasis_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SelectedONVBasis_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
