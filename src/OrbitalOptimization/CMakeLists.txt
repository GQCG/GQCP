target_sources(gqcp
        PRIVATE
            AP1roGJacobiOrbitalOptimizer.cpp
            BaseOrbitalOptimizer.cpp
            DOCINewtonOrbitalOptimizer.cpp
            JacobiOrbitalOptimizer.cpp
            JacobiRotationParameters.cpp
            NewtonOrbitalOptimizer.cpp
            OrbitalOptimizationOptions.cpp
            OrbitalRotationGenerators.cpp
            QCMethodNewtonOrbitalOptimizer.cpp
        )

add_subdirectory(Localization)