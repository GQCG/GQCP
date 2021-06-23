list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GHF_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RHF_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/UHF_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
