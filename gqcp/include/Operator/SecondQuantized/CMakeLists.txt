target_sources(gqcp
    PRIVATE
        SQHamiltonian.hpp
        SQOneElectronOperator.hpp
        SQTwoElectronOperator.hpp
        USQOneElectronOperator.hpp
        USQHamiltonian.hpp
)

add_subdirectory(ModelHamiltonian)
