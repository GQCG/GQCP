list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/DavidsonSolver_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DenseSolver_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Eigenpair_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/NewtonMinimizer_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/NewtonSystemOfEquationsSolver_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SparseSolver_test.cpp
        )

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
