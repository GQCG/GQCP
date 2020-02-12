target_sources(gqcp
    PRIVATE
        SQHamiltonian.hpp
        SQOneElectronOperator.hpp
        SQTwoElectronOperator.hpp
        USQHamiltonian.hpp
)

add_subdirectory(ModelHamiltonian)
