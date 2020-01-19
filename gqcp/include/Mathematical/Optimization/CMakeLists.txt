target_sources(gqcp
    PRIVATE
        AcceleratedIterativeSolver.hpp
        IterativeSolver.hpp
        ProxyAcceleratedIterativeSolver.hpp
)

add_subdirectory(Accelerators)
add_subdirectory(NonLinear)
add_subdirectory(Eigenproblem)
