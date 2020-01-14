list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/elements_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Molecule_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/NuclearFramework_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Nucleus_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
