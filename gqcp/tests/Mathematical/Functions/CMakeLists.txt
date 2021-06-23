list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/BoysFunction_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CartesianExponents_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CartesianGTO_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Function_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LondonCartesianGTO_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
