target_sources(gqcp
    PRIVATE
        QCMethodProtocol.hpp
        QCObjective.hpp
        QCStructure.hpp
)

add_subdirectory(Applications)
add_subdirectory(CI)
add_subdirectory(Geminals)
add_subdirectory(HF)
add_subdirectory(OrbitalOptimization)
add_subdirectory(RMP2)
