list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/AP1roG_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGBivariationalSolver_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGGeminalCoefficients_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGPSESolver_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/AP1roGVariables_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/APIGGeminalCoefficients_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/OO_AP1roG_test.cpp
        )
set(test_target_sources ${test_target_sources} PARENT_SCOPE)