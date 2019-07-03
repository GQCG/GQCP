list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/AOBasis_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CartesianExponents_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CartesianGTO_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/GTOShell_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/LibintInterfacer_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ShellSet_test.cpp
        )

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
