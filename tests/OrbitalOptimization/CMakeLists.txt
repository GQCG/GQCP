list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGJacobiOrbitalOptimizer_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DOCINewtonOrbitalOptimizer_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/JacobiRotationParameters_test.cpp
        )

add_subdirectory(Localization)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
