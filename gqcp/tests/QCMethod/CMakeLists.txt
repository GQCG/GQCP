add_subdirectory(Applications)
add_subdirectory(CI)
add_subdirectory(Geminals)
add_subdirectory(OrbitalOptimization)
add_subdirectory(RHF)
add_subdirectory(RMP2)
add_subdirectory(Utilities)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
