# In this CMake file, all CMake variables will be set


# Parse the project name into uppercase and lowercase.
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)  # Uppercase is needed in version.hpp.in
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)



# The name of the library should be equal to the project name
if(NOT LIBRARY_NAME)
    set(LIBRARY_NAME ${PROJECT_NAME})
endif()

# We want to make a static library
set(LIBRARY_TYPE STATIC)
set(EXPORT_TYPE ARCHIVE)



# Find the source folder
set(PROJECT_SOURCE_FOLDER ${CMAKE_SOURCE_DIR}/src)

# Find the source files
set(PROJECT_SOURCE_FILES
        ${PROJECT_SOURCE_FOLDER}/FockSpace/BaseFockSpace.cpp
        ${PROJECT_SOURCE_FOLDER}/FockSpace/FockSpace.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianParameters/BaseHamiltonianParameters.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianParameters/HamiltonianParameters.cpp
        ${PROJECT_SOURCE_FOLDER}/HamiltonianParameters/HamiltonianParameters_constructors.cpp
        ${PROJECT_SOURCE_FOLDER}/Operator/BaseOperator.cpp
        ${PROJECT_SOURCE_FOLDER}/Operator/OneElectronOperator.cpp
        ${PROJECT_SOURCE_FOLDER}/Operator/TwoElectronOperator.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/DIISRHFSCFSolver.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/PlainRHFSCFSolver.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/RHF.cpp
        ${PROJECT_SOURCE_FOLDER}/RHF/RHFSCFSolver.cpp
        ${PROJECT_SOURCE_FOLDER}/AOBasis.cpp
        ${PROJECT_SOURCE_FOLDER}/Atom.cpp
        ${PROJECT_SOURCE_FOLDER}/elements.cpp
        ${PROJECT_SOURCE_FOLDER}/JacobiRotationParameters.cpp
        ${PROJECT_SOURCE_FOLDER}/LibintCommunicator.cpp
        ${PROJECT_SOURCE_FOLDER}/miscellaneous.cpp
        ${PROJECT_SOURCE_FOLDER}/Molecule.cpp
        ${PROJECT_SOURCE_FOLDER}/ONV.cpp)

# Find the header folder
set(PROJECT_INCLUDE_FOLDER ${CMAKE_SOURCE_DIR}/include)

# Find the header files
set(PROJECT_INCLUDE_FILES
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/BaseFockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/FockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianParameters/BaseHamiltonianParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianParameters/HamiltonianParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianParameters/HamiltonianParameters_constructors.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/BaseOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/OneElectronOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/TwoElectronOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/DIISRHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/PlainRHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/RHF.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/RHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/AOBasis.hpp
        ${PROJECT_INCLUDE_FOLDER}/Atom.hpp
        ${PROJECT_INCLUDE_FOLDER}/common.hpp
        ${PROJECT_INCLUDE_FOLDER}/elements.hpp
        ${PROJECT_INCLUDE_FOLDER}/JacobiRotationParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/LibintCommunicator.hpp
        ${PROJECT_INCLUDE_FOLDER}/miscellaneous.hpp
        ${PROJECT_INCLUDE_FOLDER}/Molecule.hpp
        ${PROJECT_INCLUDE_FOLDER}/ONV.hpp
        ${PROJECT_INCLUDE_FOLDER}/units.hpp)

# Find the tests folder
set(PROJECT_TESTS_FOLDER ${CMAKE_SOURCE_DIR}/tests)

# Find the source files for the tests
set(PROJECT_TEST_SOURCE_FILES
        ${PROJECT_TESTS_FOLDER}/FockSpace/FockSpace_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianParameters/HamiltonianParameters_test.cpp
        ${PROJECT_TESTS_FOLDER}/HamiltonianParameters/HamiltonianParameters_constructors_test.cpp
        ${PROJECT_TESTS_FOLDER}/Operator/OneElectronOperator_test.cpp
        ${PROJECT_TESTS_FOLDER}/Operator/TwoElectronOperator_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/DIISRHFSCFSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/PlainRHFSCFSolver_test.cpp
        ${PROJECT_TESTS_FOLDER}/RHF/RHF_test.cpp
        ${PROJECT_TESTS_FOLDER}/AOBasis_test.cpp
        ${PROJECT_TESTS_FOLDER}/Atom_test.cpp
        ${PROJECT_TESTS_FOLDER}/elements_test.cpp
        ${PROJECT_TESTS_FOLDER}/JacobiRotationParameters_test.cpp
        ${PROJECT_TESTS_FOLDER}/LibintCommunicator_test.cpp
        ${PROJECT_TESTS_FOLDER}/miscellaneous_test.cpp
        ${PROJECT_TESTS_FOLDER}/Molecule_test.cpp
        ${PROJECT_TESTS_FOLDER}/ONV_test.cpp
        ${PROJECT_TESTS_FOLDER}/RMP2_test.cpp
        ${PROJECT_TESTS_FOLDER}/units_test.cpp)


# Give the user the option to specify an installation prefix. If not given as -DINSTALLATION_PREFIX, defaults to /usr/local.
if(NOT INSTALLATION_PREFIX)
    set(INSTALLATION_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()
set(PROJECT_INSTALL_DIR ${INSTALLATION_PREFIX}/${PROJECT_NAME_LOWERCASE})
set(INCLUDE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/include)
set(CMAKE_INSTALL_DIR ${PROJECT_INSTALL_DIR}/cmake)
set(LIBRARY_INSTALL_DIR ${PROJECT_INSTALL_DIR}/lib)


# Include the function that configures the executables
include(${CMAKE_SOURCE_DIR}/cmake/ConfigureExecutable.cmake)
