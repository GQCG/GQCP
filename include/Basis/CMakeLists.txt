target_sources(gqcp
        PRIVATE
            CartesianExponents.hpp
            CartesianGTO.hpp
            GTOBasisSet.hpp
            GTOShell.hpp
            ScalarBasis.hpp
            ShellSet.hpp
            SingleParticleBasis.hpp
            TransformationMatrix.hpp
        )

add_subdirectory(Integrals)
