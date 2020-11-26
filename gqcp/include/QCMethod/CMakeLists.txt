target_sources(gqcp
    PRIVATE
        QCStructure.hpp
)

add_subdirectory(CI)
add_subdirectory(Geminals)
add_subdirectory(HF)
add_subdirectory(OrbitalOptimization)
add_subdirectory(RMP2)
