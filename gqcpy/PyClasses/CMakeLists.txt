list(APPEND python_bindings_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Nucleus.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SpinorBasis.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SQOneElectronOperator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SQTwoElectronOperator.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
