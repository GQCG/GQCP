target_sources(gqcp
        PRIVATE
            AOBasis.hpp
            CartesianExponents.hpp
            CartesianGTO.hpp
            GTOBasisSet.hpp
            GTOShell.hpp
            ShellSet.hpp
        )

add_subdirectory(Integrals)
