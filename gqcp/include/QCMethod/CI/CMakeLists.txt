target_sources(gqcp
    PRIVATE
        CI.hpp
        CISolver.hpp
        DOCINewtonOrbitalOptimizer.hpp
)

add_subdirectory(HamiltonianBuilder)
