// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#pragma once


#include "Basis/Integrals/Interfaces/LibcintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintInterfacer.hpp"
#include "Basis/Integrals/Interfaces/LibintOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/Interfaces/LibintTwoElectronIntegralEngine.hpp"

#include "Basis/Integrals/BaseOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/Integrals/IntegralEngine.hpp"

#include "Basis/ScalarBasis/CartesianExponents.hpp"
#include "Basis/ScalarBasis/CartesianGTO.hpp"
#include "Basis/ScalarBasis/GTOBasisSet.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/ShellSet.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/JacobiRotationParameters.hpp"
#include "Basis/SpinorBasis/OrbitalRotationGenerators.hpp"
#include "Basis/SpinorBasis/RSpinorBasis.hpp"
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Basis/SpinorBasis/SpinComponent.hpp"
#include "Basis/SpinorBasis/USpinorBasis.hpp"

#include "Basis/transform.hpp"
#include "Basis/TransformationMatrix.hpp"

#include "FockSpace/WaveFunction/SpinUnresolvedWaveFunction.hpp"
#include "FockSpace/WaveFunction/WaveFunction.hpp"
#include "FockSpace/WaveFunction/WaveFunctionReader.hpp"

#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/BaseFrozenCoreFockSpace.hpp"
#include "FockSpace/Configuration.hpp"
#include "FockSpace/EvaluationIterator.hpp"
#include "FockSpace/FockPermutator.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceType.hpp"
#include "FockSpace/FrozenFockSpace.hpp"
#include "FockSpace/FrozenProductFockSpace.hpp"
#include "FockSpace/ONV.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"

#include "Mathematical/Algorithm/ConvergenceCriterion.hpp"
#include "Mathematical/Algorithm/IterationCycle.hpp"
#include "Mathematical/Algorithm/IterationStep.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"

#include "Mathematical/Optimization/Eigenproblem/BaseEigenproblemSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/BaseMatrixSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/DavidsonSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/DenseSolver.hpp"
#include "Mathematical/Optimization/Eigenproblem/Eigenpair.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolverOptions.hpp"
#include "Mathematical/Optimization/Eigenproblem/SparseSolver.hpp"

#include "Mathematical/Optimization/NonLinear/BaseHessianModifier.hpp"
#include "Mathematical/Optimization/NonLinear/BaseMinimizer.hpp"
#include "Mathematical/Optimization/NonLinear/IterativeIdentitiesHessianModifier.hpp"
#include "Mathematical/Optimization/NonLinear/NewtonMinimizer.hpp"
#include "Mathematical/Optimization/NonLinear/NewtonNLSystemOfEquationsSolver.hpp"
#include "Mathematical/Optimization/NonLinear/step.hpp"
#include "Mathematical/Optimization/NonLinear/UnalteringHessianModifier.hpp"

#include "Mathematical/Representation/BlockMatrix.hpp"
#include "Mathematical/Representation/BlockRankFourTensor.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "Mathematical/Representation/QCRankFourTensor.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Mathematical/Representation/SquareRankFourTensor.hpp"
#include "Mathematical/Representation/Tensor.hpp"

#include "Mathematical/CartesianDirection.hpp"
#include "Mathematical/LinearCombination.hpp"
#include "Mathematical/ScalarFunction.hpp"

#include "Molecule/elements.hpp"
#include "Molecule/Molecule.hpp"
#include "Molecule/NuclearFramework.hpp"
#include "Molecule/Nucleus.hpp"

#include "Operator/FirstQuantized/BaseFQOneElectronOperator.hpp"
#include "Operator/FirstQuantized/BaseFQTwoElectronOperator.hpp"
#include "Operator/FirstQuantized/BaseMultipoleOperator.hpp"
#include "Operator/FirstQuantized/BaseNuclearOperator.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"
#include "Operator/FirstQuantized/NuclearDipoleOperator.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"

#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"

#include "Processing/Properties/BaseElectricalResponseSolver.hpp"
#include "Processing/Properties/expectation_values.hpp"
#include "Processing/Properties/properties.hpp"
#include "Processing/Properties/RHFElectricalResponseSolver.hpp"

#include "Processing/RDM/BaseRDMBuilder.hpp"
#include "Processing/RDM/BaseSpinUnresolvedRDMBuilder.hpp"
#include "Processing/RDM/DOCIRDMBuilder.hpp"
#include "Processing/RDM/FCIRDMBuilder.hpp"
#include "Processing/RDM/FrozenCoreDOCIRDMBuilder.hpp"
#include "Processing/RDM/FrozenCoreFCIRDMBuilder.hpp"
#include "Processing/RDM/FrozenCoreRDMBuilder.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "Processing/RDM/RDMCalculator.hpp"
#include "Processing/RDM/RDMs.hpp"
#include "Processing/RDM/SelectedRDMBuilder.hpp"
#include "Processing/RDM/SpinUnresolvedFCIRDMBuilder.hpp"
#include "Processing/RDM/SpinUnresolvedRDMCalculator.hpp"
#include "Processing/RDM/TwoRDM.hpp"

#include "QCMethod/Applications/AtomicDecompositionParameters.hpp"
#include "QCMethod/Applications/DOCINewtonOrbitalOptimizer.hpp"
#include "QCMethod/Applications/DOCIRHF.hpp"
#include "QCMethod/Applications/FCI.hpp"
#include "QCMethod/Applications/FukuiDysonAnalysis.hpp"
#include "QCMethod/Applications/Hubbard.hpp"
#include "QCMethod/Applications/MullikenConstrainedFCI.hpp"

#include "QCMethod/CI/HamiltonianBuilder/DOCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreDOCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/FrozenCoreFCI.hpp"
#include "QCMethod/CI/HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "QCMethod/CI/HamiltonianBuilder/HoppingMatrix.hpp"
#include "QCMethod/CI/HamiltonianBuilder/Hubbard.hpp"
#include "QCMethod/CI/HamiltonianBuilder/SelectedCI.hpp"

#include "QCMethod/CI/CISolver.hpp"
#include "QCMethod/CI/DOCINewtonOrbitalOptimizer.hpp"

#include "QCMethod/Geminals/AP1roG.hpp"
#include "QCMethod/Geminals/AP1roGGeminalCoefficients.hpp"
#include "QCMethod/Geminals/AP1roGJacobiOrbitalOptimizer.hpp"
#include "QCMethod/Geminals/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"
#include "QCMethod/Geminals/AP1roGLagrangianOptimizer.hpp"
#include "QCMethod/Geminals/AP1roGPSEs.hpp"
#include "QCMethod/Geminals/AP1roGPSESolver.hpp"
#include "QCMethod/Geminals/APIGGeminalCoefficients.hpp"
#include "QCMethod/Geminals/GeminalCoefficientsInterface.hpp"

#include "QCMethod/OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"
#include "QCMethod/OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"

#include "QCMethod/OrbitalOptimization/BaseOrbitalOptimizer.hpp"
#include "QCMethod/OrbitalOptimization/JacobiOrbitalOptimizer.hpp"
#include "QCMethod/OrbitalOptimization/NewtonOrbitalOptimizer.hpp"
#include "QCMethod/OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"

#include "QCMethod/RHF/DIISRHFSCFSolver.hpp"
#include "QCMethod/RHF/PlainRHFSCFSolver.hpp"
#include "QCMethod/RHF/RHF.hpp"
#include "QCMethod/RHF/RHFSCFSolver.hpp"

#include "QCMethod/RMP2/RMP2.hpp"

#include "QCMethod/QCMethodProtocol.hpp"

#include "Utilities/CRTP.hpp"
#include "Utilities/linalg.hpp"
#include "Utilities/memory.hpp"
#include "Utilities/miscellaneous.hpp"
#include "Utilities/type_traits.hpp"
#include "Utilities/typedefs.hpp"
#include "Utilities/units.hpp"

#include "version.hpp"
