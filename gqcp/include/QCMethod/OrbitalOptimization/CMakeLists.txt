target_sources(gqcp
    PRIVATE
        BaseOrbitalOptimizer.hpp
        JacobiOrbitalOptimizer.hpp
        JacobiRotationParameters.hpp
        NewtonOrbitalOptimizer.hpp
        OrbitalRotationGenerators.hpp
        QCMethodNewtonOrbitalOptimizer.hpp
)

add_subdirectory(Localization)
