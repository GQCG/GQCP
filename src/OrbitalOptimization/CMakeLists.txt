target_sources(gqcp
        PRIVATE
            AP1roGJacobiOrbitalOptimizer.cpp
            BaseOrbitalOptimizer.cpp
            DOCINewtonOrbitalOptimizer.cpp
            JacobiRotationParameters.cpp
            NewtonOrbitalOptimizer.cpp
            OrbitalRotationGenerators.cpp
        )

add_subdirectory(Localization)