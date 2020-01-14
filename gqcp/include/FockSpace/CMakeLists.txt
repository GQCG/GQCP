target_sources(gqcp
    PRIVATE
        BaseFockSpace.hpp
        BaseFrozenCoreFockSpace.hpp
        Configuration.hpp
        EvaluationIterator.hpp
        FockPermutator.hpp
        FockSpace.hpp
        FockSpaceType.hpp
        FrozenFockSpace.hpp
        FrozenProductFockSpace.hpp
        ONV.hpp
        ProductFockSpace.hpp
        SelectedFockSpace.hpp
)

add_subdirectory(WaveFunction)
