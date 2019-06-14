target_sources(gqcp
        PRIVATE
            AP1roGJacobiOrbitalOptimizer.cpp
            AP1roGLagrangianNewtonOrbitalOptimizer.cpp
            BaseOrbitalOptimizer.cpp
            DOCINewtonOrbitalOptimizer.cpp
            JacobiOrbitalOptimizer.cpp
            JacobiRotationParameters.cpp
            NewtonOrbitalOptimizer.cpp
            NewtonOrbitalOptimizationOptions.cpp
            OrbitalOptimizationOptions.cpp
            OrbitalRotationGenerators.cpp
            QCMethodNewtonOrbitalOptimizer.cpp
        )

add_subdirectory(Localization)