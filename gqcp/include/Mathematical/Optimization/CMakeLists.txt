target_sources(gqcp
    PRIVATE
        ConsecutiveIteratesNormConvergence.hpp
        OptimizationEnvironment.hpp
)

add_subdirectory(Accelerator)
add_subdirectory(Eigenproblem)
add_subdirectory(Minimization)
add_subdirectory(NonLinearEquation)
