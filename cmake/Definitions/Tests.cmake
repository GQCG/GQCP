# List all library tests


# Find the tests folder
set(PROJECT_TESTS_FOLDER ${CMAKE_SOURCE_DIR}/tests)


# Find the source files for the tests
set(PROJECT_TEST_SOURCE_FILES
        ${PROJECT_TESTS_FOLDER}/CISolver/CISolver_DOCI_Davidson_test.cpp
        ${PROJECT_TESTS_FOLDER}/CISolver/CISolver_DOCI_Dense_test.cpp
        ${PROJECT_TESTS_FOLDER}/CISolver/CISolver_FCI_Davidson_test.cpp
        ${PROJECT_TESTS_FOLDER}/CISolver/CISolver_FCI_Dense_test.cpp
        ${PROJECT_TESTS_FOLDER}/CISolver/CISolver_Hubbard_Davidson_test.cpp
        ${PROJECT_TESTS_FOLDER}/CISolver/CISolver_Hubbard_Dense_test.cpp
        ${PROJECT_TESTS_FOLDER}/CISolver/CISolver_test.cpp

        ${PROJECT_TESTS_FOLDER}/FockSpace/FockSpace_test.cpp
        ${PROJECT_TESTS_FOLDER}/FockSpace/ONV_test.cpp
        ${PROJECT_TESTS_FOLDER}/FockSpace/SelectedFockSpace_test.cpp
        ${PROJECT_TESTS_FOLDER}/FockSpace/ProductFockSpace_test.cpp

        ${PROJECT_TESTS_FOLDER}/geminals/AP1roGBivariationalSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/AP1roGGeminalCoefficients_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/AP1roGPSESolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/AP1roGVariables_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/APIGGeminalCoefficients_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/OO_AP1roG_test.cpp

        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/DOCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/FCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/FrozenCoreFCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/Hubbard_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/SelectedCI_test.cpp

        ${PROJECT_TESTS_FOLDER}/HamiltonianParameters/HamiltonianParameters_test.cpp

        ${PROJECT_TESTS_FOLDER}/Localization/ERJacobiLocalizer_test.cpp
        ${PROJECT_TESTS_FOLDER}/Localization/ERNewtonLocalizer_test.cpp

        ${PROJECT_TESTS_FOLDER}/Operator/OneElectronOperator_test.cpp
        ${PROJECT_TESTS_FOLDER}/Operator/TwoElectronOperator_test.cpp

        ${PROJECT_TESTS_FOLDER}/optimization/DavidsonSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/optimization/DenseSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/optimization/Eigenpair_test.cpp
        ${PROJECT_TESTS_FOLDER}/optimization/NewtonMinimizer_test.cpp
        ${PROJECT_TESTS_FOLDER}/optimization/NewtonSystemOfEquationsSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/optimization/SparseSolver_test.cpp

        ${PROJECT_TESTS_FOLDER}/properties/expectation_values_test.cpp
        ${PROJECT_TESTS_FOLDER}/properties/properties_test.cpp

        ${PROJECT_TESTS_FOLDER}/RDM/DOCIRDMBuilder_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/FCIRDMBuilder_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/FrozenCoreRDMBuilder_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/RDMCalculator_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/SelectedRDMBuilder_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/SpinUnresolvedFCIRDMBuilder_test.cpp

        ${PROJECT_TESTS_FOLDER}/RHF/constrained_RHF_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/DIISRHFSCFSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/PlainRHFSCFSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/RHF_test.cpp

        ${PROJECT_TESTS_FOLDER}/utilities/io_test.cpp
        ${PROJECT_TESTS_FOLDER}/utilities/linalg_test.cpp
        ${PROJECT_TESTS_FOLDER}/utilities/miscellaneous_test.cpp

        ${PROJECT_TESTS_FOLDER}/WaveFunction/WaveFunction_test.cpp

        ${PROJECT_TESTS_FOLDER}/AOBasis_test.cpp
        ${PROJECT_TESTS_FOLDER}/Atom_test.cpp
        ${PROJECT_TESTS_FOLDER}/elements_test.cpp
        ${PROJECT_TESTS_FOLDER}/HoppingMatrix_test.cpp
        ${PROJECT_TESTS_FOLDER}/JacobiRotationParameters_test.cpp
        ${PROJECT_TESTS_FOLDER}/LibintCommunicator_test.cpp
        ${PROJECT_TESTS_FOLDER}/Molecule_test.cpp
        ${PROJECT_TESTS_FOLDER}/OO_DOCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/RMP2_test.cpp
        ${PROJECT_TESTS_FOLDER}/units_test.cpp
    )
