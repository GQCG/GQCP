target_sources(gqcp
        PRIVATE
            elements.cpp
            HoppingMatrix.cpp
            RMP2.cpp
        )

add_subdirectory(Basis)
add_subdirectory(CISolver)
add_subdirectory(Drivers)
add_subdirectory(FockSpace)
add_subdirectory(Geminals)
add_subdirectory(HamiltonianBuilder)
add_subdirectory(HamiltonianParameters)
add_subdirectory(Mathematical)
add_subdirectory(Molecule)
add_subdirectory(OrbitalOptimization)
add_subdirectory(Properties)
add_subdirectory(RDM)
add_subdirectory(RHF)
add_subdirectory(Utilities)
add_subdirectory(WaveFunction)

if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(Bindings)
endif()