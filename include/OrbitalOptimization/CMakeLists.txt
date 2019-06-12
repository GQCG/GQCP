target_sources(gqcp
        PRIVATE
            AP1roGJacobiOrbitalOptimizer.hpp
            BaseOrbitalOptimizer.hpp
            DOCINewtonOrbitalOptimizer.hpp
            JacobiOrbitalOptimizer.hpp
            JacobiRotationParameters.hpp
            NewtonOrbitalOptimizer.hpp
            OrbitalOptimizationOptions.hpp
            OrbitalRotationGenerators.hpp
            QCMethodNewtonOrbitalOptimizer.hpp
        )

add_subdirectory(Localization)
