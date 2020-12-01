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
#include <pybind11/pybind11.h>

namespace py = pybind11;


/**
 *  As stated in the Pybind11 FAQ (https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-reduce-the-build-time), it is good practice to split the binding of the code over multiple files. Here, we're declaring all the binding methods that we have implemented in various `_bindings.cpp` files.
 */
namespace gqcpy {


// Basis - MullikenPartitioning
void bindRMullikenPartitioning(py::module& module);
void bindUMullikenPartitioning(py::module& module);
void bindUMullikenPartitioningComponent(py::module& module);


// Basis - ScalarBasis
void bindGTOShell(py::module& module);


// Basis - SpinorBasis
void bindGSpinorBasis(py::module& module);
void bindOccupationType(py::module& module);
void bindOrbitalSpace(py::module& module);
void bindRSpinOrbitalBasis(py::module& module);
void bindUSpinOrbitalBasis(py::module& module);


// Basis - Transformations
void bindGTransformation(py::module& module);
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


// Mathematical - Algorithm
void bindAlgorithms(py::module& module);
void bindFunctionalSteps(py::module& module);
void bindIterativeAlgorithms(py::module& module);


// Mathematical - Grid
void bindCubicGrid(py::module& module);
void bindField(py::module& module);
void bindWeightedGrid(py::module& module);


// Mathematical - Optimization - Eigenproblem
void bindEigenproblemEnvironment(py::module& module);
void bindEigenproblemSolver(py::module& module);


// Mathematical - Optimization - LinearEquation
void bindLinearEquationSolver(py::module& module);


// Mathematical - Optimization - NonLinearEquation
void bindNonLinearEquationEnvironment(py::module& module);
void bindNonLinearEquationSolver(py::module& module);


// Molecule
void bindMolecule(py::module& module);
void bindNucleus(py::module& module);


// ONVBasis
void bindONVPaths(py::module& module);
void bindSeniorityZeroONVBasis(py::module& module);
void bindSpinResolvedONV(py::module& module);
void bindSpinResolvedONVBasis(py::module& module);
void bindSpinUnresolvedONV(py::module& module);
void bindSpinUnresolvedONVBasis(py::module& module);


// Operator - FirstQuantized
void bindOperator(py::module& module);


// Operator - SecondQuantized - ModelHamiltonian
void bindHoppingMatrix(py::module& module);
void bindHubbardHamiltonian(py::module& module);


// Operator - SecondQuantized
void bindGSQOneElectronOperator(py::module& module);
void bindGSQTwoElectronOperator(py::module& module);
void bindMixedUSQTwoElectronOperatorComponent(py::module& module);
void bindPureUSQTwoElectronOperatorComponent(py::module& module);
void bindRSQOneElectronOperator(py::module& module);
void bindRSQTwoElectronOperator(py::module& module);
void bindSQHamiltonians(py::module& module);
void bindUSQOneElectronOperator(py::module& module);
void bindUSQOneElectronOperatorComponent(py::module& module);
void bindUSQTwoElectronOperator(py::module& module);


// Processing - Properties
void bindDOCIElectricalResponseSolver(py::module& module);
void bindDysonOrbital(py::module& module);
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
void bindQCMethodCIs(py::module& module);


// QCMethod - Geminals
void bindAP1roGLagrangianNewtonOrbitalOptimizer(py::module& module);
void bindPSEnvironment(py::module& module);
void bindQCMethodAP1roG(py::module& module);
void bindQCMethodvAP1roG(py::module& module);


// QCMethod - HF - GHF
void bindQCMethodGHF(py::module& module);
void bindGHFSCFEnvironment(py::module& module);
void bindGHFSCFSolver(py::module& module);


// QCMethod - HF - RHF
void bindDiagonalRHFFockMatrixObjective(py::module& module);
void bindQCMethodRHF(py::module& module);
void bindRHFSCFEnvironment(py::module& module);
void bindRHFSCFSolver(py::module& module);


// QCMethod - HF - UHF
void bindQCMethodUHF(py::module& module);
void bindUHFSCFEnvironment(py::module& module);
void bindUHFSCFSolver(py::module& module);


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
void bindQCModelUHF(py::module& module);


// QuantumChemical
void bindSpin(py::module& module);


// Single includes
void bindVersion(py::module& module);


}  // namespace gqcpy


/**
 *  The actual Python binding into the gqcpy Python module.
 */
PYBIND11_MODULE(gqcpy, module) {

    // Basis - MullikenPartitioning
    gqcpy::bindRMullikenPartitioning(module);
    gqcpy::bindUMullikenPartitioning(module);
    gqcpy::bindUMullikenPartitioningComponent(module);


    // Basis - ScalarBasis
    gqcpy::bindGTOShell(module);


    // Basis - SpinorBasis
    gqcpy::bindGSpinorBasis(module);
    gqcpy::bindOccupationType(module);
    gqcpy::bindOrbitalSpace(module);
    gqcpy::bindRSpinOrbitalBasis(module);
    gqcpy::bindUSpinOrbitalBasis(module);


    // Basis - Transformations
    gqcpy::bindGTransformation(module);
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


    // Mathematical - Algorithm
    gqcpy::bindAlgorithms(module);
    gqcpy::bindFunctionalSteps(module);
    gqcpy::bindIterativeAlgorithms(module);


    // Mathematical - Grid
    gqcpy::bindCubicGrid(module);
    gqcpy::bindField(module);
    gqcpy::bindWeightedGrid(module);


    // Mathematical - Optimization - Eigenproblem
    gqcpy::bindEigenproblemEnvironment(module);
    gqcpy::bindEigenproblemSolver(module);


    // Mathematical - Optimization - LinearEquation
    gqcpy::bindLinearEquationSolver(module);


    // Mathematical - Optimization - NonLinearEquation
    gqcpy::bindNonLinearEquationEnvironment(module);
    gqcpy::bindNonLinearEquationSolver(module);


    // Molecule
    gqcpy::bindMolecule(module);
    gqcpy::bindNucleus(module);


    // ONVBasis
    gqcpy::bindONVPaths(module);
    gqcpy::bindSeniorityZeroONVBasis(module);
    gqcpy::bindSpinResolvedONV(module);
    gqcpy::bindSpinResolvedONVBasis(module);
    gqcpy::bindSpinUnresolvedONV(module);
    gqcpy::bindSpinUnresolvedONVBasis(module);


    // Operator - FirstQuantized
    gqcpy::bindOperator(module);


    // Operator - SecondQuantized - ModelHamiltonian
    gqcpy::bindHoppingMatrix(module);
    gqcpy::bindHubbardHamiltonian(module);


    // Operator - SecondQuantized
    gqcpy::bindGSQOneElectronOperator(module);
    gqcpy::bindGSQTwoElectronOperator(module);
    gqcpy::bindMixedUSQTwoElectronOperatorComponent(module);
    gqcpy::bindPureUSQTwoElectronOperatorComponent(module);
    gqcpy::bindRSQOneElectronOperator(module);
    gqcpy::bindRSQTwoElectronOperator(module);
    gqcpy::bindSQHamiltonians(module);
    gqcpy::bindUSQOneElectronOperator(module);
    gqcpy::bindUSQOneElectronOperatorComponent(module);
    gqcpy::bindUSQTwoElectronOperator(module);


    // Processing - Properties
    // gqcpy::bindDOCIElectricalResponseSolver(module);
    gqcpy::bindDysonOrbital(module);
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
    gqcpy::bindQCMethodCIs(module);


    // QCMethod - Geminals
    gqcpy::bindAP1roGLagrangianNewtonOrbitalOptimizer(module);
    gqcpy::bindPSEnvironment(module);
    gqcpy::bindQCMethodAP1roG(module);
    gqcpy::bindQCMethodvAP1roG(module);


    // QCMethod - HF - GHF
    gqcpy::bindQCMethodGHF(module);
    gqcpy::bindGHFSCFEnvironment(module);
    gqcpy::bindGHFSCFSolver(module);


    // QCMethod - HF - RHF
    gqcpy::bindDiagonalRHFFockMatrixObjective(module);
    gqcpy::bindQCMethodRHF(module);
    gqcpy::bindRHFSCFEnvironment(module);
    gqcpy::bindRHFSCFSolver(module);


    // QCMethod - HF - UHF
    gqcpy::bindQCMethodUHF(module);
    gqcpy::bindUHFSCFEnvironment(module);
    gqcpy::bindUHFSCFSolver(module);


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
    gqcpy::bindQCModelRHF(module);
    gqcpy::bindQCModelUHF(module);


    // QuantumChemical
    gqcpy::bindSpin(module);


    // Single includes
    gqcpy::bindVersion(module);
}
