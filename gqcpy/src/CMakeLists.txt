add_subdirectory(Basis)
add_subdirectory(Molecule)
add_subdirectory(Operator)
add_subdirectory(QCMethod)

set(python_test_target_sources ${python_test_target_sources} PARENT_SCOPE)
