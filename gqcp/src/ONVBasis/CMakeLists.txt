target_sources(gqcp
    PRIVATE
        BaseONVBasis.cpp
        BaseFrozenCoreONVBasis.cpp
        ONVBasis.cpp
        FrozenONVBasis.cpp
        FrozenProductONVBasis.cpp
        ONV.cpp
        ProductONVBasis.cpp
        SelectedONVBasis.cpp
)

add_subdirectory(LinearExpansion)
