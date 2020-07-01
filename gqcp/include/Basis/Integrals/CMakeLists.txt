target_sources(gqcp
    PRIVATE
        BaseOneElectronIntegralBuffer.hpp
        BaseOneElectronIntegralEngine.hpp
        BaseTwoElectronIntegralBuffer.hpp
        BaseTwoElectronIntegralEngine.hpp
        IntegralCalculator.hpp
        IntegralEngine.hpp
        McMurchieDavidsonCoefficient.hpp
        OneElectronIntegralBuffer.hpp
        OneElectronIntegralEngine.hpp
        PrimitiveDipoleIntegralEngine.hpp
        PrimitiveKineticEnergyIntegralEngine.hpp
        PrimitiveOverlapIntegralEngine.hpp
)

add_subdirectory(Interfaces)
