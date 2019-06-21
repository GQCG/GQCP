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
#ifndef GQCP_HPP
#define GQCP_HPP

#include "Basis/AOBasis.hpp"
#include "Basis/CartesianDirection.hpp"
#include "Basis/CartesianExponents.hpp"
#include "Basis/CartesianGTO.hpp"
#include "Basis/LibcintInterfacer.hpp"
#include "Basis/LibintInterfacer.hpp"
#include "Basis/Shell.hpp"
#include "Basis/ShellSet.hpp"

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
#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include "math/optimization/BaseEigenproblemSolver.hpp"
#include "math/optimization/BaseHessianModifier.hpp"
#include "math/optimization/BaseMatrixSolver.hpp"
#include "math/optimization/BaseMinimizer.hpp"
#include "math/optimization/BaseSystemOfEquationsSolver.hpp"
#include "math/optimization/DavidsonSolver.hpp"
#include "math/optimization/DenseSolver.hpp"
#include "math/optimization/Eigenpair.hpp"
#include "math/optimization/EigenproblemSolverOptions.hpp"
#include "math/optimization/IterativeIdentitiesHessianModifier.hpp"
#include "math/optimization/NewtonMinimizer.hpp"
#include "math/optimization/NewtonSystemOfEquationsSolver.hpp"
#include "math/optimization/SparseSolver.hpp"
#include "math/optimization/step.hpp"
#include "math/optimization/UnalteringHessianModifier.hpp"

#include "math/ChemicalMatrix.hpp"
#include "math/ChemicalRankFourTensor.hpp"
#include "math/LinearCombination.hpp"
#include "math/Matrix.hpp"
#include "math/ScalarFunction.hpp"
#include "math/SquareMatrix.hpp"
#include "math/SquareRankFourTensor.hpp"
#include "math/Tensor.hpp"

#include "Operator/OneElectronOperator.hpp"
#include "Operator/Operator.hpp"
#include "Operator/TwoElectronOperator.hpp"

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

#include "properties/expectation_values.hpp"
#include "properties/properties.hpp"

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

#include "utilities/linalg.hpp"
#include "utilities/miscellaneous.hpp"

#include "WaveFunction/SpinUnresolvedWaveFunction.hpp"
#include "WaveFunction/WaveFunction.hpp"
#include "WaveFunction/WaveFunctionReader.hpp"


// Single files, not in a special include directory
#include "Atom.hpp"
#include "elements.hpp"
#include "HoppingMatrix.hpp"
#include "Molecule.hpp"
#include "RMP2.hpp"
#include "typedefs.hpp"
#include "units.hpp"

#include "version.hpp"


#endif  // GQCP_HPP
