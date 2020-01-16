target_sources(gqcp
    PRIVATE
        QCMethod.hpp
        QCObjective.hpp
        QCStructure.hpp
)

add_subdirectory(Applications)
add_subdirectory(CI)
add_subdirectory(Geminals)
add_subdirectory(OrbitalOptimization)
add_subdirectory(RHF)
add_subdirectory(RMP2)
