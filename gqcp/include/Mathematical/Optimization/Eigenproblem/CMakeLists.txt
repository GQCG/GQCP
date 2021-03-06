target_sources(gqcp
    PRIVATE
        DenseDiagonalization.hpp
        Eigenpair.hpp
        EigenproblemEnvironment.hpp
        EigenproblemSolver.hpp
)

add_subdirectory(Davidson)
