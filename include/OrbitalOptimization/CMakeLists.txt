target_sources(gqcp
        PRIVATE
            AP1roGJacobiOrbitalOptimizer.hpp
            DOCINewtonOrbitalOptimizer.hpp
            OrbitalOptimizationOptions.hpp
        )

add_subdirectory(Localization)
