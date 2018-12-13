# In this CMake file, all CMake variables will be set


# Parse the project name into uppercase and lowercase.
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)  # Uppercase is needed in version.hpp.in
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)


# The name of the library should be equal to the project name
if(NOT LIBRARY_NAME)
    set(LIBRARY_NAME ${PROJECT_NAME})
endif()

# We want to make a shared/dynamic library
set(LIBRARY_TYPE SHARED)

# Find the source folder
set(PROJECT_SOURCE_FOLDER ${CMAKE_SOURCE_DIR}/src)

# Find the source files
set(PROJECT_SOURCE_FILES
        ${PROJECT_SOURCE_FOLDER}/CISolver/CISolver.cpp
        ${PROJECT_SOURCE_FOLDER}/geminals/AP1roG.cpp
        ${PROJECT_SOURCE_FOLDER}/geminals/AP1roGGeminalCoefficients.cpp
        ${PROJECT_SOURCE_FOLDER}/geminals/AP1roGJacobiOrbitalOptimizer.cpp
        ${PROJECT_SOURCE_FOLDER}/geminals/AP1roGPSESolver.cpp
        ${PROJECT_SOURCE_FOLDER}/geminals/APIGGeminalCoefficients.cpp
        ${PROJECT_SOURCE_FOLDER}/geminals/BaseAPIGGeminalCoefficients.cpp
        ${PROJECT_SOURCE_FOLDER}/FockSpace/BaseFockSpace.cpp
        ${PROJECT_SOURCE_FOLDER}/FockSpace/FockSpace.cpp
        ${PROJECT_SOURCE_FOLDER}/FockSpace/ONV.cpp
        ${PROJECT_SOURCE_FOLDER}/FockSpace/ProductFockSpace.cpp
        ${PROJECT_SOURCE_FOLDER}/FockSpace/SelectedFockSpace.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianBuilder/DOCI.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianBuilder/FCI.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianBuilder/HamiltonianBuilder.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianBuilder/Hubbard.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianParameters/BaseHamiltonianParameters.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianParameters/HamiltonianParameters.cpp
        ${PROJECT_SOURCE_FOLDER}/Localization/BaseERLocalizer.cpp
        ${PROJECT_SOURCE_FOLDER}/Localization/ERJacobiLocalizer.cpp
        ${PROJECT_SOURCE_FOLDER}/Localization/ERNewtonLocalizer.cpp
        ${PROJECT_SOURCE_FOLDER}/Operator/BaseOperator.cpp
        ${PROJECT_SOURCE_FOLDER}/Operator/OneElectronOperator.cpp
        ${PROJECT_SOURCE_FOLDER}/Operator/TwoElectronOperator.cpp
        ${PROJECT_SOURCE_FOLDER}/properties/expectation_values.cpp
        ${PROJECT_SOURCE_FOLDER}/properties/properties.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/BaseRDM.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/BaseRDMBuilder.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/DOCIRDMBuilder.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/FCIRDMBuilder.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/NRDMCalculator.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/OneRDM.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/RDMCalculator.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/RDMs.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/SelectedRDMBuilder.cpp
        ${PROJECT_SOURCE_FOLDER}/RDM/TwoRDM.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/DIISRHFSCFSolver.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/PlainRHFSCFSolver.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/RHF.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/RHFSCFSolver.cpp
        ${PROJECT_SOURCE_FOLDER}/WaveFunction/WaveFunction.cpp
        ${PROJECT_SOURCE_FOLDER}/WaveFunction/WaveFunctionReader.cpp
        ${PROJECT_SOURCE_FOLDER}/AOBasis.cpp
        ${PROJECT_SOURCE_FOLDER}/Atom.cpp
        ${PROJECT_SOURCE_FOLDER}/DOCINewtonOrbitalOptimizer.cpp
        ${PROJECT_SOURCE_FOLDER}/elements.cpp
        ${PROJECT_SOURCE_FOLDER}/HoppingMatrix.cpp
        ${PROJECT_SOURCE_FOLDER}/JacobiRotationParameters.cpp
        ${PROJECT_SOURCE_FOLDER}/LibintCommunicator.cpp
        ${PROJECT_SOURCE_FOLDER}/miscellaneous.cpp
        ${PROJECT_SOURCE_FOLDER}/Molecule.cpp
        ${PROJECT_SOURCE_FOLDER}/RMP2.cpp)

# Find the header folder
set(PROJECT_INCLUDE_FOLDER ${CMAKE_SOURCE_DIR}/include)

# Find the header files
set(PROJECT_INCLUDE_FILES
        ${PROJECT_INCLUDE_FOLDER}/CISolver/CISolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roG.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGGeminalCoefficients.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGJacobiOrbitalOptimizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGPSESolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/APIGGeminalCoefficients.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/BaseAPIGGeminalCoefficients.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/BaseFockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/Configuration.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/FockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/FockSpaceType.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/ONV.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/SelectedFockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/ProductFockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/DOCI.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/FCI.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/HamiltonianBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/Hubbard.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianParameters/BaseHamiltonianParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianParameters/HamiltonianParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/Localization/BaseERLocalizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/Localization/ERJacobiLocalizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/Localization/ERNewtonLocalizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/BaseOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/OneElectronOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/TwoElectronOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/properties/expectation_values.hpp
        ${PROJECT_INCLUDE_FOLDER}/properties/properties.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/BaseRDM.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/BaseRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/DOCIRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/FCIRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/NRDMCalculator.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/OneRDM.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/SelectedRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/RDMCalculator.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/RDMs.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/TwoRDM.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/DIISRHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/PlainRHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/RHF.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/RHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/WaveFunction/WaveFunction.hpp
        ${PROJECT_INCLUDE_FOLDER}/AOBasis.hpp
        ${PROJECT_INCLUDE_FOLDER}/Atom.hpp
        ${PROJECT_INCLUDE_FOLDER}/common.hpp
        ${PROJECT_INCLUDE_FOLDER}/DOCINewtonOrbitalOptimizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/elements.hpp
        ${PROJECT_INCLUDE_FOLDER}/HoppingMatrix.hpp
        ${PROJECT_INCLUDE_FOLDER}/JacobiRotationParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/LibintCommunicator.hpp
        ${PROJECT_INCLUDE_FOLDER}/miscellaneous.hpp
        ${PROJECT_INCLUDE_FOLDER}/Molecule.hpp
        ${PROJECT_INCLUDE_FOLDER}/RMP2.hpp
        ${PROJECT_INCLUDE_FOLDER}/units.hpp)

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
        ${PROJECT_TESTS_FOLDER}/geminals/AP1roG_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/AP1roGGeminalCoefficients_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/OO_AP1roG_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/AP1roGPSESolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/geminals/APIGGeminalCoefficients_test.cpp
        ${PROJECT_TESTS_FOLDER}/FockSpace/FockSpace_test.cpp
        ${PROJECT_TESTS_FOLDER}/FockSpace/ONV_test.cpp
        ${PROJECT_TESTS_FOLDER}/FockSpace/SelectedFockSpace_test.cpp
        ${PROJECT_TESTS_FOLDER}/FockSpace/ProductFockSpace_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/DOCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/FCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianBuilder/Hubbard_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianParameters/HamiltonianParameters_test.cpp
        ${PROJECT_TESTS_FOLDER}/Localization/ERJacobiLocalizer_test.cpp
        ${PROJECT_TESTS_FOLDER}/Localization/ERNewtonLocalizer_test.cpp
        ${PROJECT_TESTS_FOLDER}/Operator/OneElectronOperator_test.cpp
        ${PROJECT_TESTS_FOLDER}/Operator/TwoElectronOperator_test.cpp
        ${PROJECT_TESTS_FOLDER}/properties/expectation_values_test.cpp
        ${PROJECT_TESTS_FOLDER}/properties/properties_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/DOCIRDMBuilder_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/FCIRDMBuilder_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/NRDMCalculator_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/RDMCalculator_test.cpp
        ${PROJECT_TESTS_FOLDER}/RDM/SelectedRDMBuilder_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/constrained_RHF_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/DIISRHFSCFSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/PlainRHFSCFSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/RHF_test.cpp
        ${PROJECT_TESTS_FOLDER}/AOBasis_test.cpp
        ${PROJECT_TESTS_FOLDER}/Atom_test.cpp
        ${PROJECT_TESTS_FOLDER}/elements_test.cpp
        ${PROJECT_TESTS_FOLDER}/HoppingMatrix_test.cpp
        ${PROJECT_TESTS_FOLDER}/JacobiRotationParameters_test.cpp
        ${PROJECT_TESTS_FOLDER}/LibintCommunicator_test.cpp
        ${PROJECT_TESTS_FOLDER}/miscellaneous_test.cpp
        ${PROJECT_TESTS_FOLDER}/Molecule_test.cpp
        ${PROJECT_TESTS_FOLDER}/OO_DOCI_test.cpp
        ${PROJECT_TESTS_FOLDER}/RMP2_test.cpp
        ${PROJECT_TESTS_FOLDER}/units_test.cpp)

# Find the executables folder
set(PROJECT_EXECUTABLES_FOLDER ${CMAKE_SOURCE_DIR}/exe)

# Find the source files for the executables
set(PROJECT_EXE_SOURCE_FILES
    ${PROJECT_EXECUTABLES_FOLDER}/fci_lowdin.cpp
    ${PROJECT_EXECUTABLES_FOLDER}/hubbard.cpp
    ${PROJECT_EXECUTABLES_FOLDER}/oo_doci.cpp)

# Find the benchmarks folder
set(PROJECT_BENCHMARKS_FOLDER ${CMAKE_SOURCE_DIR}/benchmarks)

# Find the source files for the benchmarks
set(PROJECT_BENCH_SOURCE_FILES
        ${PROJECT_BENCHMARKS_FOLDER}/DOCI/doci_case.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/DOCI/doci_matrix.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/DOCI/doci_matvec.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/FCI/fci_hchain.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/FCI/fci_matrix.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/FCI/fci_matvec.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/Hubbard/hubbard_diagonalization.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/Hubbard/hubbard_matrix.cpp
        ${PROJECT_BENCHMARKS_FOLDER}/Hubbard/hubbard_matvec.cpp)


# Give the user the option to specify an installation prefix. If not given as -DINSTALLATION_PREFIX, defaults to /usr/local.
if(NOT INSTALLATION_PREFIX)
    set(INSTALLATION_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
set(PROJECT_INSTALL_DIR ${INSTALLATION_PREFIX}/${PROJECT_NAME_LOWERCASE})
set(INCLUDE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/include)
set(CMAKE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/cmake)
set(LIBRARY_INSTALL_DIR ${PROJECT_INSTALL_DIR}/lib)
set(BIN_INSTALL_DIR ${PROJECT_INSTALL_DIR}/bin)


# Include the function that configures target executables: tests, drivers and benchmarks
include(${CMAKE_SOURCE_DIR}/cmake/Configure/ConfigureExecutable.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/Configure/ConfigureTest.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/Configure/ConfigureDriver.cmake)
include(${CMAKE_SOURCE_DIR}/cmake/Configure/ConfigureBenchmarks.cmake)
