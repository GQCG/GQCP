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
void bindSpinorBasis(py::module& module);


// Mathematical - Algorithm
void bindAlgorithm(py::module& module);
void bindIterativeAlgorithm(py::module& module);


// Mathematical - Optimization - Eigenproblem
void bindEigenproblemEnvironment(py::module& module);
void bindEigenproblemSolver(py::module& module);


// Molecule
void bindMolecule(py::module& module);
void bindNucleus(py::module& module);


// ONVBasis
void bindSpinResolvedONVBasis(py::module& module);


// Operator
void bindSQHamiltonian(py::module& module);
void bindSQOneElectronOperator(py::module& module);
void bindSQTwoElectronOperator(py::module& module);


// QCMethod - Applications
void bindQCMethodDOCINewtonOrbitalOptimizer(py::module& module);
void bindQCMethodDOCIRHF(py::module& module);
void bindQCMethodHubbard(py::module& module);
void bindQCMethodFCI(py::module& module);
void bindQCMethodFukuiDysonAnalysis(py::module& module);
void bindMullikenConstrainedFCI(py::module& module);


// QCMethod - CI
void bindCIEnvironment(py::module& module);
void bindQCMethodCI(py::module& module);


// QCMethod - HF
void bindDiagonalRHFFockMatrixObjective(py::module& module);
void bindQCMethodRHF(py::module& module);
void bindRHFSCFEnvironment(py::module& module);
void bindRHFSCFSolver(py::module& module);


// QCModel - CI
void bindLinearExpansion(py::module& module);


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
    gqcpy::bindSpinorBasis(module);


    // Mathematical - Algorithm
    gqcpy::bindAlgorithm(module);
    gqcpy::bindIterativeAlgorithm(module);


    // Mathematical - Optimization - Eigenproblem
    gqcpy::bindEigenproblemEnvironment(module);
    gqcpy::bindEigenproblemSolver(module);


    // Molecule
    gqcpy::bindMolecule(module);
    gqcpy::bindNucleus(module);


    // ONVBasis
    gqcpy::bindSpinResolvedONVBasis(module);


    // Operator
    gqcpy::bindSQHamiltonian(module);
    gqcpy::bindSQOneElectronOperator(module);
    gqcpy::bindSQTwoElectronOperator(module);


    // QCMethod - Applications
    gqcpy::bindQCMethodDOCINewtonOrbitalOptimizer(module);
    gqcpy::bindQCMethodDOCIRHF(module);
    gqcpy::bindQCMethodHubbard(module);
    gqcpy::bindQCMethodFCI(module);
    gqcpy::bindQCMethodFukuiDysonAnalysis(module);
    gqcpy::bindMullikenConstrainedFCI(module);


    // QCMethod - CI
    gqcpy::bindCIEnvironment(module);
    gqcpy::bindQCMethodCI(module);


    // QCMethod - HF
    gqcpy::bindDiagonalRHFFockMatrixObjective(module);
    gqcpy::bindQCMethodRHF(module);
    gqcpy::bindRHFSCFEnvironment(module);
    gqcpy::bindRHFSCFSolver(module);


    // QCModel - CI
    gqcpy::bindLinearExpansion(module);


    // QCModel - HF
    gqcpy::bindQCModelRHF(module);


    // Single includes
    gqcpy::bindVersion(module);
}
