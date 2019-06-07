target_sources(gqcp
        PRIVATE
            AP1roGJacobiOrbitalOptimizer.hpp
            BaseOrbitalOptimizer.hpp
            DOCINewtonOrbitalOptimizer.hpp
            NewtonOrbitalOptimizer.hpp
            OrbitalOptimizationOptions.hpp
        )

add_subdirectory(Localization)
