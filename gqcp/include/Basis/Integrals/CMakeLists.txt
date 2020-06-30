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
        PrimitiveOverlapIntegralEngine.hpp
)

add_subdirectory(Interfaces)
