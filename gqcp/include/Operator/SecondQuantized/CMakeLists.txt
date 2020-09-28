target_sources(gqcp
    PRIVATE
        GSQOneElectronOperator.hpp
        RSQOneElectronOperator.hpp
        SimpleSQOneElectronOperator.hpp
        SQHamiltonian.hpp
        SQOneElectronOperator.hpp
        SQTwoElectronOperator.hpp
        USQOneElectronOperator.hpp
        USQTwoElectronOperator.hpp
        USQHamiltonian.hpp
)

add_subdirectory(ModelHamiltonian)
