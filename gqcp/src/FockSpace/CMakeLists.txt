target_sources(gqcp
    PRIVATE
        BaseFockSpace.cpp
        BaseFrozenCoreFockSpace.cpp
        FockSpace.cpp
        FrozenFockSpace.cpp
        FrozenProductFockSpace.cpp
        ONV.cpp
        ProductFockSpace.cpp
        SelectedFockSpace.cpp
)

add_subdirectory(WaveFunction)
