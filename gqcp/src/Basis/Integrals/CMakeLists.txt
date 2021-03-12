target_sources(gqcp
    PRIVATE
        BoysFunction.cpp
        HermiteCoulombIntegral.cpp
        IntegralEngine.cpp
        McMurchieDavidsonCoefficient.cpp
        PrimitiveAngularMomentumIntegralEngine.cpp
        PrimitiveCartesianOperatorIntegralEngine.cpp
        PrimitiveDipoleIntegralEngine.cpp
        PrimitiveKineticEnergyIntegralEngine.cpp
        PrimitiveLinearMomentumIntegralEngine.cpp
        PrimitiveOverlapIntegralEngine.cpp
)

add_subdirectory(Interfaces)
