target_sources(gqcp
    PRIVATE
        BaseOneElectronIntegralBuffer.hpp
        BaseOneElectronIntegralEngine.hpp
        BaseTwoElectronIntegralBuffer.hpp
        BaseTwoElectronIntegralEngine.hpp
        IntegralCalculator.hpp
        IntegralEngine.hpp
)

add_subdirectory(Interfaces)
