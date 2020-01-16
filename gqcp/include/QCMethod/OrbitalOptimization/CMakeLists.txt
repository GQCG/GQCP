target_sources(gqcp
    PRIVATE
        BaseOrbitalOptimizer.hpp
        JacobiOrbitalOptimizer.hpp
        NewtonOrbitalOptimizer.hpp
        QCMethodNewtonOrbitalOptimizer.hpp
)

add_subdirectory(Localization)
