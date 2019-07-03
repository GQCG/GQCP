target_sources(gqcp
        PRIVATE
            BaseOperator.hpp
            OneElectronOperator.hpp
            TwoElectronOperator.hpp
        )

add_subdirectory(FirstQuantized)
