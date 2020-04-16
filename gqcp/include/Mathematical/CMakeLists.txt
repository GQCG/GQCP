target_sources(gqcp
    PRIVATE
        CartesianDirection.hpp
        LinearCombination.hpp
        ScalarFunction.hpp
)

add_subdirectory(Algorithm)
add_subdirectory(Grid)
add_subdirectory(Optimization)
add_subdirectory(Representation)
