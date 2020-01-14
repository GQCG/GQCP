target_sources(gqcp
    PRIVATE
        gqcp.hpp
)

add_subdirectory(Basis)
add_subdirectory(FockSpace)
add_subdirectory(Mathematical)
add_subdirectory(Molecule)
add_subdirectory(Operator)
add_subdirectory(Processing)
add_subdirectory(QCMethod)
add_subdirectory(Utilities)
