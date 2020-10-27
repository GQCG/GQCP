target_sources(gqcp
    PRIVATE
        EvaluatableScalarRSQOneElectronOperator.hpp
        GSQOneElectronOperator.hpp
        RSQOneElectronOperator.hpp
        SimpleSQOneElectronOperator.hpp
        SQHamiltonian.hpp
        SQOneElectronOperator.hpp
        SQOperatorStorage.hpp
        SQOperatorStorageBase.hpp
        SQTwoElectronOperator.hpp
        USQOneElectronOperator.hpp
        USQOneElectronOperatorComponent.hpp
        USQTwoElectronOperator.hpp
)

add_subdirectory(ModelHamiltonian)
