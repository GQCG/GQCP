# List all library headers (.hpp)


# Find the header folder
set(PROJECT_INCLUDE_FOLDER ${CMAKE_SOURCE_DIR}/include)


# Find the header files
set(PROJECT_INCLUDE_FILES
        ${PROJECT_INCLUDE_FOLDER}/CISolver/CISolver.hpp

        ${PROJECT_INCLUDE_FOLDER}/FockSpace/BaseFockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/Configuration.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/FockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/FockSpaceType.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/ONV.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/SelectedFockSpace.hpp
        ${PROJECT_INCLUDE_FOLDER}/FockSpace/ProductFockSpace.hpp

        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roG.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGBivariationalSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGGeminalCoefficients.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGJacobiOrbitalOptimizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGPSESolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/AP1roGVariables.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/APIGGeminalCoefficients.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/BaseAP1roGSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/BaseAPIGVariables.hpp
        ${PROJECT_INCLUDE_FOLDER}/geminals/GeminalCoefficientsInterface.hpp

        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/DOCI.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/FCI.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/HamiltonianBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/Hubbard.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianBuilder/SelectedCI.hpp

        ${PROJECT_INCLUDE_FOLDER}/HamiltonianParameters/BaseHamiltonianParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/HamiltonianParameters/HamiltonianParameters.hpp

        ${PROJECT_INCLUDE_FOLDER}/Localization/BaseERLocalizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/Localization/ERJacobiLocalizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/Localization/ERNewtonLocalizer.hpp

        ${PROJECT_INCLUDE_FOLDER}/Operator/BaseOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/OneElectronOperator.hpp
        ${PROJECT_INCLUDE_FOLDER}/Operator/TwoElectronOperator.hpp

        ${PROJECT_INCLUDE_FOLDER}/optimization/BaseEigenproblemSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/BaseMatrixSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/BaseMinimizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/BaseSystemOfEquationsSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/DavidsonSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/DenseSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/Eigenpair.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/EigenproblemSolverOptions.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/NewtonMinimizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/NewtonSystemOfEquationsSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/SparseSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/optimization/step.hpp

        ${PROJECT_INCLUDE_FOLDER}/properties/expectation_values.hpp
        ${PROJECT_INCLUDE_FOLDER}/properties/properties.hpp

        ${PROJECT_INCLUDE_FOLDER}/RDM/BaseRDM.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/BaseRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/BaseSpinUnresolvedRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/DOCIRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/FCIRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/OneRDM.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/RDMCalculator.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/RDMs.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/SelectedRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/SpinUnresolvedFCIRDMBuilder.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/SpinUnresolvedRDMCalculator.hpp
        ${PROJECT_INCLUDE_FOLDER}/RDM/TwoRDM.hpp

        ${PROJECT_INCLUDE_FOLDER}/RHF/DIISRHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/PlainRHFSCFSolver.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/RHF.hpp
        ${PROJECT_INCLUDE_FOLDER}/RHF/RHFSCFSolver.hpp

        ${PROJECT_INCLUDE_FOLDER}/utilities/io.hpp
        ${PROJECT_INCLUDE_FOLDER}/utilities/linalg.hpp
        ${PROJECT_INCLUDE_FOLDER}/utilities/miscellaneous.hpp

        ${PROJECT_INCLUDE_FOLDER}/WaveFunction/SpinUnresolvedWaveFunction.hpp
        ${PROJECT_INCLUDE_FOLDER}/WaveFunction/WaveFunction.hpp

        ${PROJECT_INCLUDE_FOLDER}/AOBasis.hpp
        ${PROJECT_INCLUDE_FOLDER}/Atom.hpp
        ${PROJECT_INCLUDE_FOLDER}/common.hpp
        ${PROJECT_INCLUDE_FOLDER}/DOCINewtonOrbitalOptimizer.hpp
        ${PROJECT_INCLUDE_FOLDER}/elements.hpp
        ${PROJECT_INCLUDE_FOLDER}/HoppingMatrix.hpp
        ${PROJECT_INCLUDE_FOLDER}/JacobiRotationParameters.hpp
        ${PROJECT_INCLUDE_FOLDER}/LibintCommunicator.hpp
        ${PROJECT_INCLUDE_FOLDER}/Molecule.hpp
        ${PROJECT_INCLUDE_FOLDER}/RMP2.hpp
        ${PROJECT_INCLUDE_FOLDER}/units.hpp
    )
