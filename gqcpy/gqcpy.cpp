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
 *  As stated in the Pybind11 FAQ (https://pybind11.readthedocs.io/en/stable/faq.html#how-can-i-reduce-the-build-time), it is good practice to split the binding of the code over multiple files.
 */
namespace gqcpy {


// Basis
void bindGSpinorBasis(py::module& module);
void bindRSpinorBasis(py::module& module);


// Mathematical - Algorithm
void bindAlgorithm(py::module& module);
void bindIterativeAlgorithms(py::module& module);


// Mathematical - Optimization - Eigenproblem
void bindEigenproblemEnvironment(py::module& module);
void bindEigenproblemSolver(py::module& module);


// Molecule
void bindMolecule(py::module& module);
void bindNucleus(py::module& module);


// ONVBasis
void bindSeniorityZeroONVBasis(py::module& module);
void bindSpinResolvedONVBasis(py::module& module);


// Operator - FirstQuantized
void bindOperator(py::module& module);


// Operator - SecondQuantized - ModelHamiltonian
void bindHoppingMatrix(py::module& module);
void bindHubbardHamiltonian(py::module& module);


// Operator - SecondQuantized
void bindSQHamiltonian(py::module& module);
void bindSQOneElectronOperators(py::module& module);
void bindSQTwoElectronOperator(py::module& module);


// QCMethod - CI
void bindCIEnvironment(py::module& module);
void bindCIFactory(py::module& module);
void bindQCMethodCIs(py::module& module);


// QCMethod - HF
void bindDiagonalRHFFockMatrixObjective(py::module& module);
void bindQCMethodRHF(py::module& module);
void bindRHFSCFEnvironment(py::module& module);
void bindRHFSCFSolver(py::module& module);


// QCMethod
void bindQCStructures(py::module& module);


// QCModel - CI
void bindLinearExpansions(py::module& module);


// QCModel - HF
void bindQCModelRHF(py::module& module);


// Single includes
void bindVersion(py::module& module);


}  // namespace gqcpy



/**
 *  The actual Python binding into the gqcpy Python module
 */
PYBIND11_MODULE (gqcpy, module) {

    // Basis
    gqcpy::bindGSpinorBasis(module);
    gqcpy::bindRSpinorBasis(module);


    // Mathematical - Algorithm
    gqcpy::bindAlgorithm(module);
    gqcpy::bindIterativeAlgorithms(module);


    // Mathematical - Optimization - Eigenproblem
    gqcpy::bindEigenproblemEnvironment(module);
    gqcpy::bindEigenproblemSolver(module);


    // Molecule
    gqcpy::bindMolecule(module);
    gqcpy::bindNucleus(module);


    // ONVBasis
    gqcpy::bindSeniorityZeroONVBasis(module);
    gqcpy::bindSpinResolvedONVBasis(module);


    // Operator - FirstQuantized
    gqcpy::bindOperator(module);


    // Operator - SecondQuantized - ModelHamiltonian
    gqcpy::bindHoppingMatrix(module);
    gqcpy::bindHubbardHamiltonian(module);


    // Operator - SecondQuantized
    gqcpy::bindSQHamiltonian(module);
    gqcpy::bindSQOneElectronOperators(module);
    gqcpy::bindSQTwoElectronOperator(module);


    // QCMethod - CI
    gqcpy::bindCIEnvironment(module);
    gqcpy::bindCIFactory(module);
    gqcpy::bindQCMethodCIs(module);


    // QCMethod - HF
    gqcpy::bindDiagonalRHFFockMatrixObjective(module);
    gqcpy::bindQCMethodRHF(module);
    gqcpy::bindRHFSCFEnvironment(module);
    gqcpy::bindRHFSCFSolver(module);


    // QCMethod
    gqcpy::bindQCStructures(module);


    // QCModel - CI
    gqcpy::bindLinearExpansions(module);


    // QCModel - HF
    gqcpy::bindQCModelRHF(module);


    // Single includes
    gqcpy::bindVersion(module);
}
