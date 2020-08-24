target_sources(gqcp
    PRIVATE
        transform.hpp
        SpinResolvedTransformationMatrix.hpp
        TransformationMatrix.hpp
)

add_subdirectory(Integrals)
add_subdirectory(ScalarBasis)
add_subdirectory(SpinorBasis)
