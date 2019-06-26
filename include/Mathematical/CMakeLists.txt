target_sources(gqcp
        PRIVATE
            ChemicalMatrix.hpp
            ChemicalRankFourTensor.hpp
            LinearCombination.hpp
            Matrix.hpp
            ScalarFunction.hpp
            SquareMatrix.hpp
            SquareRankFourTensor.hpp
            Tensor.hpp
        )

add_subdirectory(Optimization)
