add_subdirectory(ChronusQ)

target_sources(gqcp
    PRIVATE
        LibcintInterfacer.cpp
        LibintInterfacer.cpp
)
