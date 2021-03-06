list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/GSpinorBasis_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OccupationType_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalSpace_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RSpinOrbitalBasis_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/USpinOrbitalBasis_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
