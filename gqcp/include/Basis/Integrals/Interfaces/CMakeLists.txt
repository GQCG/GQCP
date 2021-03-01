add_subdirectory(ChronusQ)

target_sources(gqcp
    PRIVATE
        LibcintInterfacer.hpp
        LibcintOneElectronIntegralBuffer.hpp
        LibcintOneElectronIntegralEngine.hpp
        LibcintTwoElectronIntegralBuffer.hpp
        LibcintTwoElectronIntegralEngine.hpp
        LibintInterfacer.hpp
        LibintOneElectronIntegralBuffer.hpp
        LibintOneElectronIntegralEngine.hpp
        LibintTwoElectronIntegralBuffer.hpp
        LibintTwoElectronIntegralEngine.hpp
)
