target_sources(gqcp
    PRIVATE
        IntegralEngine.cpp
        McMurchieDavidsonCoefficient.cpp
        PrimitiveCartesianOperatorIntegralEngine.cpp
        PrimitiveDipoleIntegralEngine.cpp
        PrimitiveKineticEnergyIntegralEngine.cpp
        PrimitiveLinearMomentumIntegralEngine.cpp
        PrimitiveOverlapIntegralEngine.cpp
)

add_subdirectory(Interfaces)
