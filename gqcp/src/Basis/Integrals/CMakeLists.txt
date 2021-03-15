target_sources(gqcp
    PRIVATE
        BoysFunction.cpp
        HermiteCoulombIntegral.cpp
        IntegralEngine.cpp
        McMurchieDavidsonCoefficient.cpp
        PrimitiveAngularMomentumIntegralEngine.cpp
        PrimitiveCartesianOperatorIntegralEngine.cpp
        PrimitiveDipoleIntegralEngine.cpp
        PrimitiveLinearMomentumIntegralEngine.cpp
)

add_subdirectory(Interfaces)
