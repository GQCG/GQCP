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


#include "CISolver/CISolver.hpp"

#include "FockSpace/BaseFockSpace.hpp"
#include "FockSpace/Configuration.hpp"
#include "FockSpace/FockSpace.hpp"
#include "FockSpace/FockSpaceType.hpp"
#include "FockSpace/ONV.hpp"
#include "FockSpace/ProductFockSpace.hpp"
#include "FockSpace/SelectedFockSpace.hpp"

#include "geminals/AP1roG.hpp"
#include "geminals/AP1roGBivariationalSolver.hpp"
#include "geminals/AP1roGGeminalCoefficients.hpp"
#include "geminals/AP1roGJacobiOrbitalOptimizer.hpp"
#include "geminals/AP1roGPSESolver.hpp"
#include "geminals/AP1roGVariables.hpp"
#include "geminals/APIGGeminalCoefficients.hpp"
#include "geminals/BaseAP1roGSolver.hpp"
#include "geminals/BaseAPIGVariables.hpp"
#include "geminals/BivariationalCoefficients.hpp"
#include "geminals/GeminalCoefficientsInterface.hpp"

#include "HamiltonianBuilder/DOCI.hpp"
#include "HamiltonianBuilder/FCI.hpp"
#include "HamiltonianBuilder/HamiltonianBuilder.hpp"
#include "HamiltonianBuilder/Hubbard.hpp"

#include "HamiltonianParameters/BaseHamiltonianParameters.hpp"
#include "HamiltonianParameters/HamiltonianParameters.hpp"

#include "Localization/BaseERLocalizer.hpp"
#include "Localization/ERJacobiLocalizer.hpp"
#include "Localization/ERNewtonLocalizer.hpp"

#include "Operator/BaseOperator.hpp"
#include "Operator/OneElectronOperator.hpp"
#include "Operator/TwoElectronOperator.hpp"

#include "optimization/BaseEigenproblemSolver.hpp"
#include "optimization/BaseMatrixSolver.hpp"
#include "optimization/BaseMinimizer.hpp"
#include "optimization/BaseSystemOfEquationsSolver.hpp"
#include "optimization/DavidsonSolver.hpp"
#include "optimization/Eigenpair.hpp"
#include "optimization/EigenproblemSolverOptions.hpp"
#include "optimization/NewtonMinimizer.hpp"
#include "optimization/NewtonSystemOfEquationsSolver.hpp"
#include "optimization/SparseSolver.hpp"
#include "optimization/step.hpp"

#include "properties/expectation_values.hpp"
#include "properties/properties.hpp"

#include "RDM/BaseRDM.hpp"
#include "RDM/BaseRDMBuilder.hpp"
#include "RDM/BaseSpinUnresolvedRDMBuilder.hpp"
#include "RDM/DOCIRDMBuilder.hpp"
#include "RDM/FCIRDMBuilder.hpp"
#include "RDM/OneRDM.hpp"
#include "RDM/RDMCalculator.hpp"
#include "RDM/RDMs.hpp"
#include "RDM/SelectedRDMBuilder.hpp"
#include "RDM/SpinUnresolvedFCIRDMBuilder.hpp"
#include "RDM/SpinUnresolvedRDMCalculator.hpp"

#include "RHF/DIISRHFSCFSolver.hpp"
#include "RHF/PlainRHFSCFSolver.hpp"
#include "RHF/RHF.hpp"
#include "RHF/RHFSCFSolver.hpp"

#include "utilities/io.hpp"
#include "utilities/linalg.hpp"
#include "utilities/miscellaneous.hpp"

#include "WaveFunction/SpinUnresolvedWaveFunction.hpp"
#include "WaveFunction/WaveFunction.hpp"
#include "WaveFunction/WaveFunctionReader.hpp"


// Single files, not in a special include directory
#include "AOBasis.hpp"
#include "Atom.hpp"
#include "typedefs.hpp"
#include "DOCINewtonOrbitalOptimizer.hpp"
#include "elements.hpp"
#include "HoppingMatrix.hpp"
#include "JacobiRotationParameters.hpp"
#include "LibintCommunicator.hpp"
#include "Molecule.hpp"
#include "OrbitalOptimizationOptions.hpp"
#include "RMP2.hpp"
#include "units.hpp"

#include "version.hpp"


#endif  // GQCP_HPP
