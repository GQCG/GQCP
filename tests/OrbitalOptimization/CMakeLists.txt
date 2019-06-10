list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/JacobiRotationParameters_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/OO_AP1roG_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/OO_DOCI_test.cpp
        )

add_subdirectory(Localization)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
