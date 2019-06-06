target_sources(gqcp
        PRIVATE
            AP1roGJacobiOrbitalOptimizer.cpp
            BaseOrbitalOptimizer.cpp
            DOCINewtonOrbitalOptimizer.cpp
        )

add_subdirectory(Localization)