target_sources(gqcp
        PRIVATE
            CartesianExponents.hpp
            CartesianGTO.hpp
            GTOBasisSet.hpp
            GTOShell.hpp
            RSpinorBasis.hpp
            ScalarBasis.hpp
            ShellSet.hpp
            TransformationMatrix.hpp
            transform.hpp
        )

add_subdirectory(Integrals)
