target_sources(gqcp
    PRIVATE
        BaseEigenproblemSolver.hpp
        BaseMatrixSolver.hpp
        DavidsonSolver.hpp
        DenseDiagonalization.hpp
        DenseSolver.hpp
        Eigenpair.hpp
        EigenproblemEnvironment.hpp
        EigenproblemSolver.hpp
        EigenproblemSolverOptions.hpp
        SparseSolver.hpp
)

add_subdirectory(Davidson)
