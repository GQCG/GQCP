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

#include "Basis/CartesianExponents.hpp"
#include "Basis/CartesianGTO.hpp"
#include "Basis/GTOBasisSet.hpp"
#include "Basis/GTOShell.hpp"
#include "Basis/ScalarBasis.hpp"
#include "Basis/ShellSet.hpp"
#include "Basis/SingleParticleBasis.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Basis/transform.hpp"

#include "Basis/Integrals/BaseOneElectronIntegralBuffer.hpp"
#include "Basis/Integrals/BaseOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralBuffer.hpp"
#include "Basis/Integrals/BaseTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/Integrals/IntegralEngine.hpp"

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

#include "CISolver/CISolver.hpp"

#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/BaseFrozenCoreFockSpace.hpp"
#include "FockSpace/Configuration.hpp"
#include "FockSpace/EvaluationMatrix.hpp"
#include "FockSpace/FockPermutator.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceType.hpp"
#include "FockSpace/FrozenFockSpace.hpp"
#include "FockSpace/FrozenProductFockSpace.hpp"
#include "FockSpace/ONV.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"

#include "Geminals/AP1roG.hpp"
#include "Geminals/AP1roGGeminalCoefficients.hpp"
#include "Geminals/AP1roGLagrangianOptimizer.hpp"
#include "Geminals/AP1roGPSEs.hpp"
#include "Geminals/AP1roGPSESolver.hpp"
#include "Geminals/AP1roGVariables.hpp"
#include "Geminals/APIGGeminalCoefficients.hpp"
#include "Geminals/BaseAP1roGSolver.hpp"
#include "Geminals/BaseAPIGVariables.hpp"
#include "Geminals/GeminalCoefficientsInterface.hpp"

#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/FrozenCoreCI.hpp"
#include "HamiltonianBuilder/FrozenCoreDOCI.hpp"
#include "HamiltonianBuilder/FrozenCoreFCI.hpp"
#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"
#include "HamiltonianBuilder/SelectedCI.hpp"

#include "HamiltonianParameters/AtomicDecompositionParameters.hpp"

#include "Mathematical/Optimization/BaseEigenproblemSolver.hpp"
#include "Mathematical/Optimization/BaseHessianModifier.hpp"
#include "Mathematical/Optimization/BaseMatrixSolver.hpp"
#include "Mathematical/Optimization/BaseMinimizer.hpp"
#include "Mathematical/Optimization/DavidsonSolver.hpp"
#include "Mathematical/Optimization/DenseSolver.hpp"
#include "Mathematical/Optimization/Eigenpair.hpp"
#include "Mathematical/Optimization/EigenproblemSolverOptions.hpp"
#include "Mathematical/Optimization/IterativeIdentitiesHessianModifier.hpp"
#include "Mathematical/Optimization/NewtonMinimizer.hpp"
#include "Mathematical/Optimization/NewtonNLSystemOfEquationsSolver.hpp"
#include "Mathematical/Optimization/SparseSolver.hpp"
#include "Mathematical/Optimization/step.hpp"
#include "Mathematical/Optimization/UnalteringHessianModifier.hpp"

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

#include "OrbitalOptimization/Localization/ERJacobiLocalizer.hpp"
#include "OrbitalOptimization/Localization/ERNewtonLocalizer.hpp"

#include "OrbitalOptimization/AP1roGJacobiOrbitalOptimizer.hpp"
#include "OrbitalOptimization/AP1roGLagrangianNewtonOrbitalOptimizer.hpp"
#include "OrbitalOptimization/BaseOrbitalOptimizer.hpp"
#include "OrbitalOptimization/DOCINewtonOrbitalOptimizer.hpp"
#include "OrbitalOptimization/JacobiOrbitalOptimizer.hpp"
#include "OrbitalOptimization/JacobiRotationParameters.hpp"
#include "OrbitalOptimization/NewtonOrbitalOptimizer.hpp"
#include "OrbitalOptimization/OrbitalRotationGenerators.hpp"
#include "OrbitalOptimization/QCMethodNewtonOrbitalOptimizer.hpp"

#include "Properties/expectation_values.hpp"
#include "Properties/properties.hpp"

#include "QCMethod/DOCINewtonOrbitalOptimizer.hpp"
#include "QCMethod/FCI.hpp"
#include "QCMethod/Hubbard.hpp"
#include "QCMethod/MullikenConstrainedFCI.hpp"

#include "RDM/BaseRDMBuilder.hpp"
#include "RDM/BaseSpinUnresolvedRDMBuilder.hpp"
#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"
#include "RDM/FrozenCoreDOCIRDMBuilder.hpp"
#include "RDM/FrozenCoreFCIRDMBuilder.hpp"
#include "RDM/FrozenCoreRDMBuilder.hpp"
#include "RDM/OneRDM.hpp"
#include "RDM/RDMCalculator.hpp"
#include "RDM/RDMs.hpp"
#include "RDM/SelectedRDMBuilder.hpp"
#include "RDM/SpinUnresolvedFCIRDMBuilder.hpp"
#include "RDM/SpinUnresolvedRDMCalculator.hpp"
#include "RDM/TwoRDM.hpp"

#include "RHF/DIISRHFSCFSolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "RHF/RHF.hpp"
#include "RHF/RHFSCFSolver.hpp"

#include "Utilities/linalg.hpp"
#include "Utilities/miscellaneous.hpp"

#include "WaveFunction/SpinUnresolvedWaveFunction.hpp"
#include "WaveFunction/WaveFunction.hpp"
#include "WaveFunction/WaveFunctionReader.hpp"


// Single files, not in a special include directory
#include "elements.hpp"
#include "HoppingMatrix.hpp"
#include "RMP2.hpp"
#include "typedefs.hpp"
#include "units.hpp"

#include "version.hpp"
