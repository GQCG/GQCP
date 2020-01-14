target_sources(gqcp
    PRIVATE
        BaseOrbitalOptimizer.cpp
        JacobiOrbitalOptimizer.cpp
        NewtonOrbitalOptimizer.cpp
        QCMethodNewtonOrbitalOptimizer.cpp
)

add_subdirectory(Localization)
