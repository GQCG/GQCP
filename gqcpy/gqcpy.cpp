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
#include "Mathematical/Representation/Tensor.hpp"

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


/**
 *  As stated in the Pybind11 FAQ (https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-reduce-the-build-time), it is good practice to split the binding of the code over multiple files. Here, we're declaring all the binding methods that we have implemented in various `_bindings.cpp` files.
 */
namespace gqcpy {


// Basis - NonOrthogonalBasis
void bindGLowdinPairingBases(py::module& module);
void bindRLowdinPairingBases(py::module& module);
void bindULowdinPairingBases(py::module& module);


// Basis - Integrals
void bindFunctionalOneElectronPrimitiveIntegralEngine(py::module& module);
void bindFunctionalOneElectronIntegralEngine(py::module& module);
void bindFunctionalTwoElectronPrimitiveIntegralEngine(py::module& module);
void bindFunctionalTwoElectronIntegralEngine(py::module& module);
void bindHermiteCoulombIntegral(py::module& module);
void bindMcMurchieDavidsonCoefficient(py::module& module);
void bindIntegralEngine(py::module& module);


// Basis - NonOrthogonalBasis
void bindGNonOrthogonalStateBases(py::module& module);
void bindRNonOrthogonalStateBases(py::module& module);
void bindUNonOrthogonalStateBases(py::module& module);


// Basis - ScalarBasis
void bindGTOShell(py::module& module);
void bindScalarBasis(py::module& module);
void bindShellSet(py::module& module);


// Basis - SpinorBasis
void bindCurrentDensityMatrixElement(py::module& module);
void bindGSpinorBases(py::module& module);
void bindOccupationType(py::module& module);
void bindOrbitalSpace(py::module& module);
void bindRSpinOrbitalBases(py::module& module);
void bindUSpinOrbitalBases(py::module& module);


// Basis - Transformations
void bindGTransformations(py::module& module);
void bindRTransformation(py::module& module);
void bindUTransformation(py::module& module);
void bindUTransformationComponent(py::module& module);


// DensityMatrix
void bindG1DM(py::module& module);
void bindG2DM(py::module& module);
void bindMixedSpinResolved2DMComponent(py::module& module);
void bindOrbital1DM(py::module& module);
void bindOrbital2DM(py::module& module);
void bindPureSpinResolved2DMComponent(py::module& module);
void bindSpinDensity1DM(py::module& module);
void bindSpinResolved1DM(py::module& module);
void bindSpinResolved1DMComponent(py::module& module);
void bindSpinResolved2DM(py::module& module);

// Domain - HubbardDomain
void bindDiscreteDomain(py::module& module);
void bindHubbardDomain(py::module& module);


// Domain - MullikenDomain
void bindGMullikenDomain(py::module& module);
void bindRMullikenDomain(py::module& module);
void bindUMullikenDomain(py::module& module);
void bindUMullikenDomainComponent(py::module& module);


// Mathematical - Algorithm
void bindAlgorithms(py::module& module);
void bindFunctionalSteps(py::module& module);
void bindIterativeAlgorithms(py::module& module);


// Mathematical - Functions
void bindCartesianDirection(py::module& module);
void bindCartesianExponents(py::module& module);
void bindCartesianGTO(py::module& module);


// Mathematical - Grid
void bindCubicGrid(py::module& module);
void bindField(py::module& module);
void bindWeightedGrid(py::module& module);


// Mathematical - Optimization - Eigenproblem
void bindEigenproblemEnvironment(py::module& module);
void bindEigenproblemSolver(py::module& module);
void bindGeneralizedEigenproblemEnvironment(py::module& module);
void bindGeneralizedEigenproblemSolver(py::module& module);


// Mathematical - Optimization - LinearEquation
void bindLinearEquationSolver(py::module& module);


// Mathematical - Optimization - NonLinearEquation
void bindNonLinearEquationEnvironment(py::module& module);
void bindNonLinearEquationSolver(py::module& module);


// Molecule
void bindMolecule(py::module& module);
void bindNuclearFramework(py::module& module);
void bindNucleus(py::module& module);


// ONVBasis
void bindONVPaths(py::module& module);
void bindSeniorityZeroONVBasis(py::module& module);
void bindSpinResolvedONV(py::module& module);
void bindSpinResolvedONVBasis(py::module& module);
void bindSpinResolvedOperatorString(py::module& module);
void bindSpinResolvedSelectedONVBasis(py::module& module);
void bindSpinUnresolvedONV(py::module& module);
void bindSpinUnresolvedOperatorString(py::module& module);
void bindSpinUnresolvedONVBasis(py::module& module);
void bindSpinUnresolvedSelectedONVBasis(py::module& module);


// Operator - FirstQuantized
void bindAngularMomentumOperator(py::module& module);
void bindCoulombRepulsionOperator(py::module& module);
void bindCurrentDensityOperator(py::module& module);
void bindDiamagneticOperator(py::module& module);
void bindElectronicDensityOperator(py::module& module);
void bindElectronicDipoleOperator(py::module& module);
void bindElectronicSpin_zOperator(py::module& module);
void bindElectronicQuadrupoleOperator(py::module& module);
void bindElectronicSpinOperator(py::module& module);
void bindElectronicSpinSquaredOperator(py::module& module);
void bindFQMolecularHamiltonian(py::module& module);
void bindFQMolecularMagneticHamiltonian(py::module& module);
void bindFQMolecularPauliHamiltonian(py::module& module);
void bindKineticOperator(py::module& module);
void bindLinearMomentumOperator(py::module& module);
void bindNuclearAttractionOperator(py::module& module);
void bindNuclearDipoleOperator(py::module& module);
void bindNuclearRepulsionOperator(py::module& module);
void bindOrbitalZeemanOperator(py::module& module);
void bindOverlapOperator(py::module& module);
void bindSpinZeemanOperator(py::module& module);


// Operator - SecondQuantized - ModelHamiltonian
void bindAdjacencyMatrix(py::module& module);
void bindHoppingMatrix(py::module& module);
void bindHubbardHamiltonian(py::module& module);


// Operator - SecondQuantized
void bindEvaluableRSQOneElectronOperator(py::module& module);
void bindGSQOneElectronOperator(py::module& module);
void bindGSQTwoElectronOperator(py::module& module);
void bindMixedUSQTwoElectronOperatorComponent(py::module& module);
void bindPureUSQTwoElectronOperatorComponent(py::module& module);
void bindRSQOneElectronOperator(py::module& module);
void bindRSQTwoElectronOperator(py::module& module);
void bindScalarGSQOneElectronOperatorProduct(py::module& module);
void bindScalarUSQOneElectronOperatorProduct(py::module& module);
void bindSQHamiltonians(py::module& module);
void bindUSQOneElectronOperator(py::module& module);
void bindUSQOneElectronOperatorComponent(py::module& module);
void bindUSQTwoElectronOperator(py::module& module);


// Partition
void bindDiscreteDomainPartition(py::module& module);
void bindONVPartition(py::module& module);
void bindSpinResolvedElectronPartition(py::module& module);
void bindSpinUnresolvedElectronPartition(py::module& module);


// Physical
void bindHomogeneousMagneticField(py::module& module);


// Processing - Properties
void bindDOCIElectricalResponseSolver(py::module& module);
void bindDysonOrbital(py::module& module);
void bindExpectationValues(py::module& module);
void bindRHFElectricalResponseSolver(py::module& module);
void bindvAP1roGElectricalResponseSolver(py::module& module);


// QCMethod - CC
void bindQCMethodCCD(py::module& module);
void bindQCMethodCCSD(py::module& module);
void bindCCSDEnvironment(py::module& module);
void bindCCDSolver(py::module& module);
void bindCCSDSolver(py::module& module);


// QCMethod - CI
void bindCIEnvironments(py::module& module);
void bindCIFactory(py::module& module);
void bindDOCINewtonOrbitalOptimizers(py::module& module);
void bindDOCINewtonOrbitalOptimizerFactory(py::module& module);
void bindQCMethodCI(py::module& module);


// QCMethod - Geminals
void bindAP1roGLagrangianNewtonOrbitalOptimizer(py::module& module);
void bindPSEnvironment(py::module& module);
void bindQCMethodAP1roG(py::module& module);
void bindQCMethodvAP1roG(py::module& module);


// QCMethod - HF - GHF
void bindQCMethodsGHF(py::module& module);
void bindGHFSCFEnvironments(py::module& module);
void bindGHFSCFSolvers(py::module& module);


// QCMethod - HF - RHF
void bindDiagonalRHFFockMatrixObjectives(py::module& module);
void bindQCMethodsRHF(py::module& module);
void bindRHFSCFEnvironments(py::module& module);
void bindRHFSCFSolvers(py::module& module);


// QCMethod - HF - UHF
void bindQCMethodsUHF(py::module& module);
void bindUHFSCFEnvironments(py::module& module);
void bindUHFSCFSolvers(py::module& module);


// QCMethod - NOCI
void bindNOCIEnvironments(py::module& module);
void bindNOCIFactory(py::module& module);
void bindQCMethodNOCI(py::module& module);


// QCMethod
void bindQCStructures(py::module& module);


// QCModel - CC
void bindQCModelCCD(py::module& module);
void bindQCModelCCSD(py::module& module);
void bindT1Amplitudes(py::module& module);
void bindT2Amplitudes(py::module& module);


// QCModel - CI
void bindLinearExpansions(py::module& module);


// QCModel - Geminals
void bindAP1roGGeminalCoefficients(py::module& module);
void bindQCModelAP1roG(py::module& module);
void bindQCModelvAP1roG(py::module& module);


// QCModel - HF
void bindQCModelGHF(py::module& module);
void bindQCModelRHF(py::module& module);
void bindQCModelsUHF(py::module& module);


// QCModel - HF - StabilityMatrices
void bindGHFStabilityMatrices(py::module& module);
void bindRHFStabilityMatrices(py::module& module);
void bindUHFStabilityMatrices(py::module& module);


// QCModel - NOCI
void bindNOCIExpansions(py::module& module);


// QuantumChemical
void bindSpin(py::module& module);
void bindSpinResolvedTypes(py::module& module);


// Single includes
void bindVersion(py::module& module);


}  // namespace gqcpy


/**
 *  The actual Python binding into the gqcpy Python module.
 */
PYBIND11_MODULE(gqcpy, module) {

    // Basis - BiOrthogonalBasis
    gqcpy::bindGLowdinPairingBases(module);
    gqcpy::bindRLowdinPairingBases(module);
    gqcpy::bindULowdinPairingBases(module);


    // Basis - Integrals
    gqcpy::bindFunctionalOneElectronPrimitiveIntegralEngine(module);
    gqcpy::bindFunctionalOneElectronIntegralEngine(module);
    gqcpy::bindFunctionalTwoElectronPrimitiveIntegralEngine(module);
    gqcpy::bindFunctionalTwoElectronIntegralEngine(module);
    gqcpy::bindHermiteCoulombIntegral(module);
    gqcpy::bindMcMurchieDavidsonCoefficient(module);
    gqcpy::bindIntegralEngine(module);


    // Basis - NonOrthogonalBasis
    gqcpy::bindGNonOrthogonalStateBases(module);
    gqcpy::bindRNonOrthogonalStateBases(module);
    gqcpy::bindUNonOrthogonalStateBases(module);


    // Basis - ScalarBasis
    gqcpy::bindGTOShell(module);
    gqcpy::bindScalarBasis(module);
    gqcpy::bindShellSet(module);


    // Basis - SpinorBasis
    gqcpy::bindCurrentDensityMatrixElement(module);
    gqcpy::bindGSpinorBases(module);
    gqcpy::bindOccupationType(module);
    gqcpy::bindOrbitalSpace(module);
    gqcpy::bindRSpinOrbitalBases(module);
    gqcpy::bindUSpinOrbitalBases(module);


    // Basis - Transformations
    gqcpy::bindGTransformations(module);
    gqcpy::bindRTransformation(module);
    gqcpy::bindUTransformation(module);
    gqcpy::bindUTransformationComponent(module);


    // DensityMatrix
    gqcpy::bindG1DM(module);
    gqcpy::bindG2DM(module);
    gqcpy::bindMixedSpinResolved2DMComponent(module);
    gqcpy::bindOrbital1DM(module);
    gqcpy::bindOrbital2DM(module);
    gqcpy::bindPureSpinResolved2DMComponent(module);
    gqcpy::bindSpinDensity1DM(module);
    gqcpy::bindSpinResolved1DM(module);
    gqcpy::bindSpinResolved1DMComponent(module);
    gqcpy::bindSpinResolved2DM(module);


    // Domain - HubbardDomain
    gqcpy::bindDiscreteDomain(module);
    gqcpy::bindHubbardDomain(module);


    // Domain - MullikenDomain
    gqcpy::bindGMullikenDomain(module);
    gqcpy::bindRMullikenDomain(module);
    gqcpy::bindUMullikenDomain(module);
    gqcpy::bindUMullikenDomainComponent(module);


    // Mathematical - Algorithm
    gqcpy::bindAlgorithms(module);
    gqcpy::bindFunctionalSteps(module);
    gqcpy::bindIterativeAlgorithms(module);


    // Mathematical - Functions
    gqcpy::bindCartesianDirection(module);
    gqcpy::bindCartesianExponents(module);
    gqcpy::bindCartesianGTO(module);


    // Mathematical - Grid
    gqcpy::bindCubicGrid(module);
    gqcpy::bindField(module);
    gqcpy::bindWeightedGrid(module);


    // Mathematical - Optimization - Eigenproblem
    gqcpy::bindEigenproblemEnvironment(module);
    gqcpy::bindEigenproblemSolver(module);
    gqcpy::bindGeneralizedEigenproblemEnvironment(module);
    gqcpy::bindGeneralizedEigenproblemSolver(module);


    // Mathematical - Optimization - LinearEquation
    gqcpy::bindLinearEquationSolver(module);


    // Mathematical - Optimization - NonLinearEquation
    gqcpy::bindNonLinearEquationEnvironment(module);
    gqcpy::bindNonLinearEquationSolver(module);


    // Molecule
    gqcpy::bindMolecule(module);
    gqcpy::bindNuclearFramework(module);
    gqcpy::bindNucleus(module);


    // ONVBasis
    gqcpy::bindONVPaths(module);
    gqcpy::bindSeniorityZeroONVBasis(module);
    gqcpy::bindSpinResolvedONV(module);
    gqcpy::bindSpinResolvedONVBasis(module);
    gqcpy::bindSpinResolvedOperatorString(module);
    gqcpy::bindSpinResolvedSelectedONVBasis(module);
    gqcpy::bindSpinUnresolvedONV(module);
    gqcpy::bindSpinUnresolvedOperatorString(module);
    gqcpy::bindSpinUnresolvedONVBasis(module);
    gqcpy::bindSpinUnresolvedSelectedONVBasis(module);


    // Operator - FirstQuantized
    gqcpy::bindAngularMomentumOperator(module);
    gqcpy::bindCoulombRepulsionOperator(module);
    gqcpy::bindCurrentDensityOperator(module);
    gqcpy::bindDiamagneticOperator(module);
    gqcpy::bindElectronicDensityOperator(module);
    gqcpy::bindElectronicDipoleOperator(module);
    gqcpy::bindElectronicQuadrupoleOperator(module);
    gqcpy::bindElectronicSpin_zOperator(module);
    gqcpy::bindElectronicSpinOperator(module);
    gqcpy::bindElectronicSpinSquaredOperator(module);
    gqcpy::bindKineticOperator(module);
    gqcpy::bindLinearMomentumOperator(module);
    gqcpy::bindFQMolecularHamiltonian(module);
    gqcpy::bindFQMolecularMagneticHamiltonian(module);
    gqcpy::bindFQMolecularPauliHamiltonian(module);
    gqcpy::bindNuclearAttractionOperator(module);
    gqcpy::bindNuclearDipoleOperator(module);
    gqcpy::bindNuclearRepulsionOperator(module);
    gqcpy::bindOrbitalZeemanOperator(module);
    gqcpy::bindOverlapOperator(module);
    gqcpy::bindSpinZeemanOperator(module);


    // Operator - SecondQuantized - ModelHamiltonian
    gqcpy::bindAdjacencyMatrix(module);
    gqcpy::bindHoppingMatrix(module);
    gqcpy::bindHubbardHamiltonian(module);


    // Operator - SecondQuantized
    gqcpy::bindEvaluableRSQOneElectronOperator(module);
    gqcpy::bindGSQOneElectronOperator(module);
    gqcpy::bindGSQTwoElectronOperator(module);
    gqcpy::bindMixedUSQTwoElectronOperatorComponent(module);
    gqcpy::bindPureUSQTwoElectronOperatorComponent(module);
    gqcpy::bindRSQOneElectronOperator(module);
    gqcpy::bindRSQTwoElectronOperator(module);
    gqcpy::bindScalarGSQOneElectronOperatorProduct(module);
    gqcpy::bindScalarUSQOneElectronOperatorProduct(module);
    gqcpy::bindSQHamiltonians(module);
    gqcpy::bindUSQOneElectronOperator(module);
    gqcpy::bindUSQOneElectronOperatorComponent(module);
    gqcpy::bindUSQTwoElectronOperator(module);


    // Partition
    gqcpy::bindDiscreteDomainPartition(module);
    gqcpy::bindONVPartition(module);
    gqcpy::bindSpinResolvedElectronPartition(module);
    gqcpy::bindSpinUnresolvedElectronPartition(module);


    // Physical
    gqcpy::bindHomogeneousMagneticField(module);


    // Processing - Properties
    // gqcpy::bindDOCIElectricalResponseSolver(module);
    gqcpy::bindDysonOrbital(module);
    gqcpy::bindExpectationValues(module);
    gqcpy::bindRHFElectricalResponseSolver(module);
    gqcpy::bindvAP1roGElectricalResponseSolver(module);


    // QCMethod - CC
    gqcpy::bindQCMethodCCD(module);
    gqcpy::bindQCMethodCCSD(module);
    gqcpy::bindCCSDEnvironment(module);
    gqcpy::bindCCDSolver(module);
    gqcpy::bindCCSDSolver(module);


    // QCMethod - CI
    gqcpy::bindCIEnvironments(module);
    gqcpy::bindCIFactory(module);
    gqcpy::bindDOCINewtonOrbitalOptimizers(module);
    gqcpy::bindDOCINewtonOrbitalOptimizerFactory(module);
    gqcpy::bindQCMethodCI(module);


    // QCMethod - Geminals
    gqcpy::bindAP1roGLagrangianNewtonOrbitalOptimizer(module);
    gqcpy::bindPSEnvironment(module);
    gqcpy::bindQCMethodAP1roG(module);
    gqcpy::bindQCMethodvAP1roG(module);


    // QCMethod - HF - GHF
    gqcpy::bindQCMethodsGHF(module);
    gqcpy::bindGHFSCFEnvironments(module);
    gqcpy::bindGHFSCFSolvers(module);


    // QCMethod - HF - RHF
    gqcpy::bindDiagonalRHFFockMatrixObjectives(module);
    gqcpy::bindQCMethodsRHF(module);
    gqcpy::bindRHFSCFEnvironments(module);
    gqcpy::bindRHFSCFSolvers(module);


    // QCMethod - HF - UHF
    gqcpy::bindQCMethodsUHF(module);
    gqcpy::bindUHFSCFEnvironments(module);
    gqcpy::bindUHFSCFSolvers(module);


    // QCMethod - NOCI
    gqcpy::bindNOCIEnvironments(module);
    gqcpy::bindNOCIFactory(module);
    gqcpy::bindQCMethodNOCI(module);


    // QCMethod
    gqcpy::bindQCStructures(module);


    // QCModel - CC
    gqcpy::bindQCModelCCD(module);
    gqcpy::bindQCModelCCSD(module);
    gqcpy::bindT1Amplitudes(module);
    gqcpy::bindT2Amplitudes(module);


    // QCModel - CI
    gqcpy::bindLinearExpansions(module);


    // QCModel - Geminals
    gqcpy::bindAP1roGGeminalCoefficients(module);
    gqcpy::bindQCModelAP1roG(module);
    gqcpy::bindQCModelvAP1roG(module);


    // QCModel - HF
    gqcpy::bindQCModelGHF(module);
    gqcpy::bindQCModelRHF(module);
    gqcpy::bindQCModelsUHF(module);


    // QCModel - HF - StabilityMatrices
    gqcpy::bindGHFStabilityMatrices(module);
    gqcpy::bindRHFStabilityMatrices(module);
    gqcpy::bindUHFStabilityMatrices(module);


    // QCModel - NOCI
    gqcpy::bindNOCIExpansions(module);


    // QuantumChemical
    gqcpy::bindSpin(module);
    gqcpy::bindSpinResolvedTypes(module);


    // Single includes
    gqcpy::bindVersion(module);
}
