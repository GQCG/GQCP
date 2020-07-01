target_sources(gqcp
    PRIVATE
        IntegralEngine.cpp
        McMurchieDavidsonCoefficient.cpp
        PrimitiveDipoleIntegralEngine.cpp
        PrimitiveKineticEnergyIntegralEngine.cpp
        PrimitiveOverlapIntegralEngine.cpp
)

add_subdirectory(Interfaces)
