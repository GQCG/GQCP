target_sources(gqcp
    PRIVATE
        BaseOrbitalOptimizer.cpp
        JacobiOrbitalOptimizer.cpp
        JacobiRotationParameters.cpp
        NewtonOrbitalOptimizer.cpp
        OrbitalRotationGenerators.cpp
        QCMethodNewtonOrbitalOptimizer.cpp
)

add_subdirectory(Localization)
