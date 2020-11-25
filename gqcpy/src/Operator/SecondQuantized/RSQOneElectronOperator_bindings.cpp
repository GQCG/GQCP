// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Register multiple variants of `RSQOneElectronOperator` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which the `RSQOneElectronOperator`s should be registered.
 */
void bindRSQOneElectronOperator(py::module& module) {

    // Define Python classes related to `RSQOneElectronOperator` and expose their interfaces.
    py::class_<ScalarRSQOneElectronOperator<double>> py_ScalarRSQOneElectronOperator_d {module, "ScalarRSQOneElectronOperator_d", "A (real) restricted one-electron operator, which is suited for expressing non-relativistic (spin-free) one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_ScalarRSQOneElectronOperator_d);
    bindScalarSQOneElectronOperatorParameterInterface(py_ScalarRSQOneElectronOperator_d);


    py::class_<ScalarRSQOneElectronOperator<complex>> py_ScalarRSQOneElectronOperator_cd {module, "ScalarRSQOneElectronOperator_cd", "A (complex) restricted one-electron operator, which is suited for expressing non-relativistic (spin-free) one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_ScalarRSQOneElectronOperator_cd);
    bindScalarSQOneElectronOperatorParameterInterface(py_ScalarRSQOneElectronOperator_cd);


    py::class_<VectorRSQOneElectronOperator<double>> py_VectorRSQOneElectronOperator_d {module, "VectorRSQOneElectronOperator_d", "A (real) restricted one-electron operator, which is suited for expressing non-relativistic (spin-free) one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_VectorRSQOneElectronOperator_d);


    py::class_<VectorRSQOneElectronOperator<complex>> py_VectorRSQOneElectronOperator_cd {module, "VectorRSQOneElectronOperator_cd", "A (complex) restricted one-electron operator, which is suited for expressing non-relativistic (spin-free) one-electron operators."};

    bindSimpleSQOneElectronOperatorInterface(py_VectorRSQOneElectronOperator_cd);
}


}  // namespace gqcpy
