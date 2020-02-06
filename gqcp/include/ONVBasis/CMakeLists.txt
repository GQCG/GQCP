target_sources(gqcp
    PRIVATE
        BaseONVBasis.hpp
        BaseFrozenCoreONVBasis.hpp
        Configuration.hpp
        EvaluationIterator.hpp
        ONVManipulator.hpp
        ONVBasis.hpp
        ONVBasisType.hpp
        FrozenONVBasis.hpp
        FrozenProductONVBasis.hpp
        ONV.hpp
        ProductONVBasis.hpp
        SelectedONVBasis.hpp
)

add_subdirectory(LinearExpansion)
