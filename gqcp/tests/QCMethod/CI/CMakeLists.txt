list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DOCI_test.cpp
    # ${CMAKE_CURRENT_SOURCE_DIR}/DOCINewtonOrbitalOptimizer_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FCI_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Hubbard_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/selected_CI_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
