list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DOCI_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FCI_builder_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FrozenCoreDOCI_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FrozenCoreFCI_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Hubbard_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SelectedCI_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
