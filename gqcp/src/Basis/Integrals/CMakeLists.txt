target_sources(gqcp
    PRIVATE
        IntegralEngine.cpp
        McMurchieDavidsonCoefficient.cpp
        PrimitiveOverlapIntegralEngine.cpp
        PrimitiveKineticEnergyIntegralEngine.cpp
)

add_subdirectory(Interfaces)
