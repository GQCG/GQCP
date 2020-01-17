target_sources(gqcp
    PRIVATE
        IterativeSolver.hpp
)

add_subdirectory(NonLinear)
add_subdirectory(Eigenproblem)
